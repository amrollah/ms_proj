function Yq = m5ppredict(model, Xq)
% m5ppredict
% Predicts response values for the given query points Xq using M5' tree or
% ensemble of trees.
%
% Call:
%   Yq = m5ppredict(model, Xq)
%
% Input:
%   model         : M5' model or a cell array of M5' models, if ensemble of
%                   trees is to be used.
%   Xq            : A matrix of query data points. Missing values in Xq
%                   must be indicated as NaN.
%
% Output:
%   Yq            : A column vector of predicted response values. If model
%                   is an ensemble, Yq is a matrix whose rows correspond to
%                   Xq rows (i.e., observations) and columns correspond to
%                   each ensemble size (i.e., the increasing number of
%                   trees), the values in the very last column being the
%                   values for a full ensemble.
%
% Remarks:
% 1. If the data contains categorical variables with more than two
%    categories, they are transformed in a number of synthetic binary
%    variables in exactly the same way as m5pbuild does it.
% 2. Every previously unseen value of a binary or categorical variable is
%    treated as NaN.

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

if nargin < 2
    error('Not enough input arguments.');
end
[nq, m] = size(Xq);

numModels = length(model);
if (numModels > 1)
    models = model;
    model = models{1};
else
    if iscell(model)
        model = model{1};
    end
end

Yq = zeros(nq,numModels);

% Transform all the categorical variables to binary ones (exactly the same way as with training data)
if any(model.binCat.binCat >= 2)
    binCatCounter = 0;
    synthCounter = 1;
    Xnew = NaN(nq, model.binCat.varMap{end}(end));
    for i = 1 : m
        if model.binCat.binCat(i) > 2
            binCatCounter = binCatCounter + 1;
            %{
            % Warn if a value was not seen when the tree was built
            XX = Xq(:,i);
            XX = unique(XX(~isnan(XX))); % no NaNs
            diff = setdiff(XX, model.binCat.catVals{binCatCounter});
            if ~isempty(diff)
                fprintf('Warning: Categorical variable x%d has one or more previously unseen values: %s. Treating as NaN.', i, mat2str(diff(:)'));
            end
            %}
            len = length(model.binCat.catVals{binCatCounter});
            for j = 1 : len
                where = Xq(:,i) == model.binCat.catVals{binCatCounter}(j);
                Xnew(where,synthCounter:(synthCounter-1 + j-1)) = 1;
                Xnew(where,(synthCounter-1 + j):(synthCounter-1 + len-1)) = 0;
            end
            synthCounter = synthCounter + len-1;
        elseif model.binCat.binCat(i) == 2
            binCatCounter = binCatCounter + 1;
            %{
            % Warn if a value was not seen when the tree was built
            XX = Xq(:,i);
            XX = unique(XX(~isnan(XX))); % no NaNs
            diff = setdiff(XX, model.binCat.catVals{binCatCounter});
            if ~isempty(diff)
                fprintf('Warning: Binary variable x%d has one or more previously unseen values: %s. Treating as NaN.', i, mat2str(diff(:)'));
            end
            %}
            val = Xq(:,i);
            catVals = model.binCat.catVals{binCatCounter};
            where = (val == catVals(1)) | (val == catVals(2));
            Xnew(where,synthCounter) = val(where,1);
            synthCounter = synthCounter + 1;
        else
            Xnew(:,synthCounter) = Xq(:,i);
            synthCounter = synthCounter + 1;
        end
    end
    Xq = Xnew;
end

if numModels == 1
    for i = 1 : nq
        Yq(i) = predictsingle(model.tree, Xq(i,:), model.trainParams.modelTree);
    end
else
    % rows correspond to observations, columns correspond to increasing number of trees
    for j = 1 : numModels
        for i = 1 : nq
            Yq(i,j:numModels) = Yq(i,j:numModels) + ...
                predictsingle(models{j}.tree, Xq(i,:), model.trainParams.modelTree);
        end
        Yq(:,j) = Yq(:,j) ./ j;
    end
end
return
