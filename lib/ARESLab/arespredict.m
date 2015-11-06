function Yq = arespredict(model, Xq)
% arespredict
% Predicts response values for the given query points Xq using ARES model.
%
% Call:
%   Yq = arespredict(model, Xq)
%
% Input:
%   model         : ARES model or a cell array of ARES models (for
%                   multi-response modelling).
%   Xq            : A matrix of query data points.
%
% Output:
%   Yq            : Either a column vector of predicted response values or,
%                   for multi-response modelling, a matrix with columns
%                   corresponding to response variables.

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

% Last update: September 7, 2015

if nargin < 2
    error('Not enough input arguments.');
end

numModels = length(model);
if (numModels > 1)
    models = model;
    model = models{1};
else
    if iscell(model)
        model = model{1};
    end
end

X = ones(size(Xq,1),length(model.knotdims)+1);
if model.trainParams.cubic
    for i = 1 : length(model.knotdims)
        X(:,i+1) = createbasisfunction(Xq, X, model.knotdims{i}, model.knotsites{i}, ...
            model.knotdirs{i}, model.parents(i), model.minX, model.maxX, model.t1(i,:), model.t2(i,:));
    end
else
    for i = 1 : length(model.knotdims)
        X(:,i+1) = createbasisfunction(Xq, X, model.knotdims{i}, model.knotsites{i}, ...
            model.knotdirs{i}, model.parents(i), model.minX, model.maxX);
    end
end

if numModels == 1
    Yq = X * model.coefs;
else
    Yq = zeros(size(Xq, 1), numModels);
    for k = 1 : numModels
        Yq(:,k) = X * models{k}.coefs;
    end
end
return
