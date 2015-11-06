function results = arestest(model, Xtst, Ytst, weights)
% arestest
% Tests ARES model on a test data set (Xtst, Ytst).
%
% Call:
%   results = arestest(model, Xtst, Ytst, weights)
%
% Input:
%   model         : ARES model or a cell array of ARES models (for
%                   multi-response modelling).
%   Xtst, Ytst    : Xtst is a matrix with rows corresponding to testing
%                   observations, and columns corresponding to input
%                   variables. Ytst is either a column vector of response
%                   values or, for multi-response data, a matrix with
%                   columns corresponding to response variables.
%   weights       : Optional. A vector of weights for observations. See
%                   description for function aresbuild.
%
% Output:
%   results       : A structure of different error measures calculated on
%                   the test data set. For multi-response data, all error
%                   measures are given for each model separately in a row
%                   vector. The structure has the following fields:
%     MAE         : Mean Absolute Error.
%     MSE         : Mean Squared Error.
%     RMSE        : Root Mean Squared Error.
%     RRMSE       : Relative Root Mean Squared Error.
%     R2          : Coefficient of Determination.

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

if nargin < 3
    error('Not enough input arguments.');
end
if isempty(Xtst) || isempty(Ytst)
    error('Data is empty.');
end
[n, dy] = size(Ytst); % number of observations and number of output variables
if (size(Xtst, 1) ~= n)
    error('The number of rows in Xtst and Ytst should be equal.');
end
if (nargin < 4)
    weights = [];
else
    if (~isempty(weights)) && ...
       ((size(weights,1) ~= n) || (size(weights,2) ~= 1))
        error('weights vector is of wrong size.');
    end
end

numModels = length(model);
if numModels == 1
    if dy ~= 1
        error('Ytst should have one column.');
    end
    residuals = arespredict(model, Xtst) - Ytst;
    if isempty(weights)
        results.MAE = mean(abs(residuals));
        results.MSE = mean(residuals .^ 2);
    else
        results.MAE = sum(abs(residuals) .* weights) / sum(weights);
        results.MSE = sum(residuals .^ 2 .* weights) / sum(weights);
    end
    results.RMSE = sqrt(results.MSE);
    if n > 1
        variance = var(Ytst, 1);
        results.RRMSE = results.RMSE / sqrt(variance);
        results.R2 = 1 - results.MSE / variance;
    else
        results.RRMSE = Inf;
        results.R2 = -Inf;
    end
else
    if dy ~= numModels
        error('The number of columns in Ytst should match the number of models.');
    end
    results.MAE = Inf(1,dy);
    results.MSE = Inf(1,dy);
    results.RMSE = Inf(1,dy);
    results.RRMSE = Inf(1,dy);
    results.R2 = -Inf(1,dy);
    Yq = arespredict(model, Xtst);
    for i = 1 : dy
        residuals = Yq(:,i) - Ytst(:,i);
        if isempty(weights)
            results.MAE(i) = mean(abs(residuals));
            results.MSE(i) = mean(residuals .^ 2);
        else
            results.MAE(i) = sum(abs(residuals) .* weights) / sum(weights);
            results.MSE(i) = sum(residuals .^ 2 .* weights) / sum(weights);
        end
        results.RMSE(i) = sqrt(results.MSE(i));
        if n > 1
            variance = var(Ytst(:,i), 1);
            results.RRMSE(i) = results.RMSE(i) / sqrt(variance);
            results.R2(i) = 1 - results.MSE(i) / variance;
        end
    end
end

return
