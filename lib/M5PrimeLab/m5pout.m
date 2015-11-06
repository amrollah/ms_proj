function m5pout(model, showNumCases, precision, plotTree, plotFontSize, dealWithNaN)
% m5pout
% Prints or plots M5' tree in a human-readable form. Does not work with
% ensembles.
%
% Call:
%   m5pout(model, showNumCases, precision, plotTree, plotFontSize, dealWithNaN)
%
% All the input arguments, except the first one, are optional. Empty values
% are also accepted (the corresponding default values will be used).
%
% Input:
%   model         : M5' model.
%   showNumCases  : Whether to show the number of training observations
%                   corresponding to each leaf (default value = true).
%   precision     : Number of digits in the model coefficients, split
%                   values etc. (default value = 15)
%   plotTree      : Whether to plot the tree instead of printing it
%                   (default value = false). In the plotted tree, left
%                   child of a node corresponds to outcome 'true' and right
%                   child to 'false'.
%   plotFontSize  : Font size for text in the plot (default value = 10).
%   dealWithNaN   : Whether to display how the tree deals with unknown
%                   values (NaN, displayed as '?'). (default value =
%                   false)
%
% Remarks:
% 1. For smoothed trees, the smoothing process is already done in m5pbuild,
%    therefore if you want to see unsmoothed versions (which are usually
%    easier to interpret) you should build trees with smoothing disabled.
% 2. If the training data has categorical variables with more than two
%    categories, the corresponding synthetic binary variables are shown.

% =========================================================================
% M5PrimeLab: M5' regression tree, model tree, and tree ensemble toolbox for Matlab/Octave
% Author: Gints Jekabsons (gints.jekabsons@rtu.lv)
% URL: http://www.cs.rtu.lv/jekabsons/
%
% Copyright (C) 2010-2015  Gints Jekabsons
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
% =========================================================================

% Last update: October 6, 2015

if nargin < 1
    error('Not enough input arguments.');
end
if length(model) > 1
    error('This function works with single trees only.');
else
    if iscell(model)
        model = model{1};
    end
end
if (nargin < 2) || isempty(showNumCases)
    showNumCases = true;
end
if (nargin < 3) || isempty(precision)
    precision = 15;
end
if (nargin < 4) || isempty(plotTree)
    plotTree = false;
end
if (nargin < 5) || isempty(plotFontSize)
    plotFontSize = 10;
end
if (nargin < 6) || isempty(dealWithNaN)
    dealWithNaN = false;
end

% Show synthetic variables
if any(model.binCat.binCat > 2)
    disp('Synthetic variables:');
    indCounter = 0;
    binCatCounter = 0;
    for i = 1 : length(model.binCat.binCat)
        if model.binCat.binCat(i) > 2
            binCatCounter = binCatCounter + 1;
            for j = 1 : length(model.binCat.catVals{binCatCounter})-1
                indCounter = indCounter + 1;
                str = num2str(model.binCat.catVals{binCatCounter}(j+1:end)', [' %.' num2str(precision) 'g,']);
                disp(['z' num2str(indCounter) ' = 1, if x' num2str(i) ' is in {' str(1:end-1) '} else = 0']);
            end
        else
            if model.binCat.binCat(i) == 2
                binCatCounter = binCatCounter + 1;
            end
            indCounter = indCounter + 1;
            disp(['z' num2str(indCounter) ' = x' num2str(i)]);
        end
    end
    zx = 'z';
else
    zx = 'x';
end

if isfield(model.binCat, 'minVals')
    minVals = model.binCat.minVals;
else
    minVals = [];
end
if ~plotTree
    if model.trainParams.smoothingK > 0
        disp('The tree (smoothed):');
    else
        disp('The tree:');
    end
    output(model.tree, model.trainParams.modelTree, model.binCat.binCatNew, ...
           minVals, 0, zx, showNumCases, precision, dealWithNaN);
    printinfo(model);
else
    if model.trainParams.modelTree
        if model.trainParams.smoothingK > 0
            disp('Models (smoothed):');
        else
            disp('Models:');
        end
    end
    figure('color', [1,1,1]);
    axis off;
    hold on;
    len = []; % number of nodes in a column
    analyzeChildren(model.tree, 1);
    pos = zeros(max(len), length(len)); % positions of nodes in columns
    for i = 1 : length(len)
        num = len(i);
        step = 1 / (num + 1);
        where = 1:-step:0;
        pos(1:num,i) = where(2:end-1);
    end
    idx = zeros(1, length(len)); % index of current positions in columns
    p = ['%.' num2str(precision) 'g'];
    numModel = 1;
    plotChildren(model.tree, 0, 1, model.trainParams.modelTree, model.binCat.binCatNew, plotFontSize);
end

    function analyzeChildren(node, depth)
        if length(len) >= depth
            len(depth) = len(depth) + 1;
        else
            len(depth) = 1;
        end
        if ~node.interior
            return;
        end
        if (~isempty(node.left))
            analyzeChildren(node.left, depth + 1);
        end
        if (~isempty(node.right))
            analyzeChildren(node.right, depth + 1);
        end
    end

    function plotChildren(node, x, depth, modelTree, binCatNew, plotFontSize)
        idx(depth) = idx(depth) + 1;
        myY = pos(idx(depth), depth);
        
        if ~node.interior
            if ~modelTree
                text(-myY, x - 1.6, num2str(node.value,p), ...
                    'FontSize', plotFontSize, 'HorizontalAlignment', 'center');
            else
                % show regression model
                str = ['M' num2str(numModel) ' = ' num2str(node.modelCoefs(1),p)];
                for k = 1 : length(node.modelAttrIdx)
                    if node.modelCoefs(k+1) >= 0
                        str = [str ' +'];
                    else
                        str = [str ' '];
                    end
                    str = [str num2str(node.modelCoefs(k+1),p) '*' zx num2str(node.modelAttrIdx(k))];
                end
                if dealWithNaN && (~isempty(node.modelAttrIdx))
                    str = [str ' (replace unknown: '];
                    for k = 1 : length(node.modelAttrIdx)
                        if k > 1
                            str = [str ', '];
                        end
                        str = [str zx num2str(node.modelAttrIdx(k)) '=' num2str(node.modelAttrAvg(k),p)];
                    end
                    str = [str ')'];
                end
                disp(str);
                str = ['M' num2str(numModel)];
                text(-myY, x - 1.65, str, 'FontSize', plotFontSize, 'HorizontalAlignment', 'center');
                numModel = numModel + 1;
            end
            if showNumCases
                text(-myY, x - 3.65, ['(' num2str(node.numCases) ')'], ...
                    'FontSize', plotFontSize, 'HorizontalAlignment', 'center');
            end
            if (depth == 1) % if the tree has only one node - the root node
                plot(-[myY;myY], [x+10;x+10], 'w');
                plot(-[myY;myY], [x;x], 'k-o');
                plot(-[myY;myY], [x-10;x-10], 'w');
            end
            return;
        end
        
        if binCatNew(node.splitAttr) % a binary variable (might be synthetic)
            str = ([zx num2str(node.splitAttr) '==' num2str(minVals(node.splitAttr),p)]);
            if dealWithNaN && node.nanLeft
                if depth == 1
                    % one line for root, so that the text does not go outside the image
                    str = [str ' OR ' zx num2str(node.splitAttr) '==?'];
                else
                    text(-myY, x + 4, str, 'FontSize', plotFontSize, 'HorizontalAlignment', 'center');
                    str = ['OR ' zx num2str(node.splitAttr) '==?'];
                end
            end
            text(-myY, x + 2, str, 'FontSize', plotFontSize, 'HorizontalAlignment', 'center');
        else % a continuous variable
            str = ([zx num2str(node.splitAttr) '<=']);
            if depth == 1
                % one line for root, so that the text does not go outside the image
                str = [str num2str(node.splitLocation)];
                if dealWithNaN && node.nanLeft
                    str = [str ' OR ' zx num2str(node.splitAttr) '==?'];
                end
            else
                if dealWithNaN && node.nanLeft
                    str = [str num2str(node.splitLocation)];
                    text(-myY, x + 4, str, 'FontSize', plotFontSize, 'HorizontalAlignment', 'center');
                    str = ['OR ' zx num2str(node.splitAttr) '==?'];
                else
                    text(-myY, x + 4, str, 'FontSize', plotFontSize, 'HorizontalAlignment', 'center');
                    str = num2str(node.splitLocation);
                end
            end
        end
        text(-myY, x + 2, str, 'FontSize', plotFontSize, 'HorizontalAlignment', 'center');
        newX = x - 10;
        if (~isempty(node.left))
            newY = pos(idx(depth + 1) + 1, depth + 1);
            plot(-[myY;newY], [x;newX], 'k-o');
            plotChildren(node.left, newX, depth + 1, modelTree, binCatNew, plotFontSize);
        end
        if (~isempty(node.right))
            newY = pos(idx(depth + 1) + 1, depth + 1);
            plot(-[myY;newY], [x;newX], 'k-o');
            plotChildren(node.right, newX, depth + 1, modelTree, binCatNew, plotFontSize);
        end
    end

end

function output(node, modelTree, binCatNew, minVals, offset, zx, showNumCases, precision, dealWithNaN)
p = ['%.' num2str(precision) 'g'];
if node.interior
    if binCatNew(node.splitAttr) % a binary variable (might be synthetic)
        str = [repmat(' ',1,offset) 'if ' zx num2str(node.splitAttr) ' == ' num2str(minVals(node.splitAttr),p)];
    else % a continuous variable
        str = [repmat(' ',1,offset) 'if ' zx num2str(node.splitAttr) ' <= ' num2str(node.splitLocation)];
    end
    if dealWithNaN && node.nanLeft
        str = [str ' OR ' zx num2str(node.splitAttr) ' == ?'];
    end
    disp(str);
    output(node.left, modelTree, binCatNew, minVals, offset + 1, zx, showNumCases, precision, dealWithNaN);
    disp([repmat(' ',1,offset) 'else']);
    output(node.right, modelTree, binCatNew, minVals, offset + 1, zx, showNumCases, precision, dealWithNaN);
    %disp([repmat(' ',1,offset) 'end']);
else
    if ~modelTree
        str = [repmat(' ',1,offset) 'y = ' num2str(node.value,p)];
    else
        if dealWithNaN
            for k = 1 : length(node.modelAttrIdx)
                disp([repmat(' ',1,offset) 'if ' zx num2str(node.modelAttrIdx(k)) ' == ?, ' zx num2str(node.modelAttrIdx(k)) ' = ' num2str(node.modelAttrAvg(k),p)]);
            end
        end
        % show regression model
        str = [repmat(' ',1,offset) 'y = ' num2str(node.modelCoefs(1),p)];
        for k = 1 : length(node.modelAttrIdx)
            if node.modelCoefs(k+1) >= 0
                str = [str ' +'];
            else
                str = [str ' '];
            end
            str = [str num2str(node.modelCoefs(k+1),p) '*' zx num2str(node.modelAttrIdx(k))];
        end
    end
    if showNumCases
        str = [str ' (' num2str(node.numCases) ')'];
    end
    disp(str);
end
end
