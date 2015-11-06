function aresplot(model, minX, maxX, vals, gridSize)
% aresplot
% Plots ARES model. For datsets with one input variable, plots the
% function together with its knot locations. For datasets with more than
% one input variable, plots the surface.
% The function works with single-response models only.
%
% Call:
%   aresplot(model, minX, maxX, vals, gridSize)
%
% All the input arguments, except the first one, are optional. Empty values
% are also accepted (the corresponding default values will be used).
%
% Input:
%   model         : ARES model.
%   minX, maxX    : User-defined minimum and maximum values for each input
%                   variable (this is the same type of data as in
%                   model.minX and model.maxX). If not supplied, the
%                   model.minX and model.maxX values will be used.
%   vals          : Only used when the number of input variables is larger
%                   than 2. This is a vector of fixed values for all the
%                   input variables except the two varied in the plot. The
%                   two varied variables are identified using NaN values.
%                   By default the two first variables will be varied and
%                   all the other will be fixed at (min+max)/2.
%   gridSize      : Grid size for the surface. The default value is 100 for
%                   datasets with one input variable and 50 for datasets
%                   with two input variables.

% =========================================================================
% ARESLab: Adaptive Regression Splines toolbox for Matlab/Octave
% Author: Gints Jekabsons (gints.jekabsons@rtu.lv)
% URL: http://www.cs.rtu.lv/jekabsons/
%
% Copyright (C) 2009-2015  Gints Jekabsons
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

% Last update: September 30, 2015

if nargin < 1
    error('Not enough input arguments.');
end
if length(model) > 1
    error('This function works with single-response models only.');
else
    if iscell(model)
        model = model{1};
    end
end
if (nargin < 2) || isempty(minX)
    minX = model.minX;
else
    if length(minX) ~= length(model.minX)
        error('Vector minX is of wrong size.');
    end
end
if (nargin < 3) || isempty(maxX)
    maxX = model.maxX;
else
    if length(maxX) ~= length(model.maxX)
        error('Vector maxX is of wrong size.');
    end
end
if (nargin < 5) || isempty(gridSize)
    if length(model.minX) < 2
        gridSize = 100;
    else
        gridSize = 50;
    end
end

if length(model.minX) < 2
    if model.trainParams.cubic
        step = (maxX - minX) / gridSize;
        X = (minX:step:maxX)';
    else
        X = [minX; maxX]; % for piecewise-linear models there are no need in more than this
    end
    knotsX = [];
    nBasis = length(model.knotsites);
    if nBasis > 0
        % Find all knot locations
        for i = 1 : nBasis
            newKnotX = model.knotsites{i};
            if (newKnotX > minX) && (newKnotX < maxX)
                knotsX = union(knotsX, newKnotX);
            end
        end
        knotsX = unique(knotsX);
        knotsX = knotsX(:);
        X = union(X, knotsX); % adding knot locations to X so that the graph shows them correctly and accurately
    end
    figure;
    Y = arespredict(model, X);
    plot(X, Y);
    grid on;
    xlabel('x');
    ylabel('y');
    if ~isempty(knotsX)
        % Add all knots to the plot
        hold on;
        knotsY = arespredict(model, knotsX);
        plot(knotsX, knotsY, 'ko', 'MarkerSize', 8);
        ylim = get(gca, 'ylim');
        for i = 1 : length(knotsX)
            plot([knotsX(i) knotsX(i)], [ylim(1) knotsY(i)], 'k--');
        end
        hold off;
    end
    return
end

if (nargin < 4) || isempty(vals)
    % By default we will plot by varying the first two variables
    ind1 = 1;
    ind2 = 2;
    vals = NaN(1,length(minX));
    if length(minX) > 2
        vals(3:end) = (minX(3:end) + maxX(3:end)) ./ 2;
    end
else
    if length(minX) ~= length(vals)
        error('Vector vals is of wrong size.');
    end
    tmp = 0;
    for i = 1 : length(vals)
        if isnan(vals(i))
            if tmp == 0
                ind1 = i;
                tmp = 1;
            elseif tmp == 1
                ind2 = i;
                tmp = 2;
            else
                tmp = 3;
                break;
            end
        end
    end
    if tmp ~= 2
        error('Vector vals should contain NaN exactly two times.');
    end
end

% Creating grid
step1 = (maxX(ind1) - minX(ind1)) / gridSize;
step2 = (maxX(ind2) - minX(ind2)) / gridSize;
[X1,X2] = meshgrid(minX(ind1):step1:maxX(ind1), minX(ind2):step2:maxX(ind2));
XX1 = reshape(X1, numel(X1), 1);
XX2 = reshape(X2, numel(X2), 1);

% Fill the other columns with the fixed values
X = zeros(size(XX1,1), length(minX));
X(:,ind1) = XX1;
X(:,ind2) = XX2;
for i = 1 : length(minX)
    if (i ~= ind1) && (i ~= ind2)
        X(:,i) = vals(i);
    end
end

% Calculate Y and plot the surface
YY = arespredict(model, X);
Y = reshape(YY, size(X1,1), size(X2,2));
figure;
surfc(X1, X2, Y);
xlabel(['x' num2str(ind1)]);
ylabel(['x' num2str(ind2)]);
zlabel('y');
return
