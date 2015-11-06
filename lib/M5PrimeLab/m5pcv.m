function [resultsTotal, resultsFolds] = m5pcv(X, Y, trainParams, ...
    isBinCat, k, shuffle, trainParamsEnsemble, verbose)
% m5pcv
% Tests M5' performance using k-fold Cross-Validation.
%
% Call:
%   [resultsTotal, resultsFolds] = ...
%       m5pcv(X, Y, trainParams, isBinCat, k, shuffle, ...
%             trainParamsEnsemble, verbose)
%
% All the input arguments, except the first two, are optional. Empty values
% are also accepted (the corresponding default values will be used).
% For more stable results, call m5pcv a few times and average the results.
% Note that, if parameter shuffle is set to true, this function employs
% random number generator for which you can set seed before calling the
% function, if you require it.
%
% Input:
%   X, Y          : Observations. Missing values in X must be indicated as
%                   NaN. (see function m5pbuild for details)
%   trainParams   : A structure of training parameters. If not provided,
%                   default values will be used (see function m5pparams for
%                   details).
%   isBinCat      : See description for function m5pbuild.
%   k             : Value of k for k-fold Cross-Validation. The typical
%                   values are 5 or 10. For Leave-One-Out Cross-Validation
%                   set k equal to n. (default value = 10)
%   shuffle       : Whether to shuffle the order of observations before
%                   performing Cross-Validation. (default value = true)
%   trainParamsEnsemble : A structure of parameters for building ensembles
%                   of trees. If not provided, a single tree is built. See
%                   function m5pparamsensemble for details.
%   verbose       : Whether to output additional information to console.
%                   (default value = true)
%
% Output:
%   resultsTotal  : A structure of results averaged over Cross-Validation
%                   folds. For tree ensembles, the structure contains
%                   fields that are column vectors with one value for each
%                   ensemble size, the very last value being value for a
%                   full ensemble.
%   resultsFolds  : A structure of row vectors of results for each Cross-
%                   Validation fold. For tree ensembles, the structure
%                   contains matrices whose rows correspond to Cross-
%                   Validation folds while columns correspond to each
%                   ensemble size, the very last value being a value for a
%                   full ensemble.
%   Both structures have the following fields:
%     MAE         : Mean Absolute Error.
%     MSE         : Mean Squared Error.
%     RMSE        : Root Mean Squared Error.
%     RRMSE       : Relative Root Mean Squared Error. Not reported for
%                   Leave-One-Out Cross-Validation.
%     R2          : Coefficient of Determination. Not reported for
%                   Leave-One-Out Cross-Validation.
%     nRules      : Number of rules in tree. For ensembles of trees this
%                   field is omitted.
%     nVars       : Number of input variables included in tree. This counts
%                   original variables (not synthetic ones that are
%                   automatically made). For ensembles of trees this field
%                   is omitted.

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

% Last update: September 26, 2015

if nargin < 2
    error('Not enough input arguments.');
end

if isempty(X) || isempty(Y)
    error('Data is empty.');
end
[n, d] = size(X); % number of observations and number of input variables
if size(Y,1) ~= n
    error('The number of rows in X and Y should be equal.');
end
if size(Y,2) ~= 1
    error('Y should have one column.');
end

if nargin < 3
    trainParams = [];
end
if nargin < 4
    isBinCat = [];
end
if (nargin < 5) || isempty(k)
    k = 10;
end
if k < 2
    error('k should not be smaller than 2.');
end
if k > n
    error('k should not be larger than the number of observations.');
end
if (nargin < 6) || isempty(shuffle)
    shuffle = true;
end
if (nargin < 7) || isempty(trainParamsEnsemble)
    trainParamsEnsemble = [];
else
    % turning off some features so that no time is wasted
    trainParamsEnsemble.getOOBError = false;
    trainParamsEnsemble.getVarImportance = 0;
    trainParamsEnsemble.verboseNumIter = 0;
end
if (nargin < 8) || isempty(verbose)
    verbose = true;
end

if shuffle
    ind = randperm(n); % shuffle the data
else
    ind = 1 : n;
end

% divide the data into k subsets (for compatibility with Octave, not using special Matlab functions)
minsize = floor(n / k);
sizes = repmat(minsize, k, 1);
remainder = n - minsize * k;
if remainder > 0
    sizes(1:remainder) = minsize + 1;
end
offsets = ones(k, 1);
for i = 2 : k
    offsets(i) = offsets(i-1) + sizes(i-1);
end

if isempty(trainParamsEnsemble)
    resultsFolds.nRules = NaN(1,k);
    resultsFolds.nVars = zeros(1,k);
    resultsFolds.MAE = Inf(1,k);
    resultsFolds.MSE = Inf(1,k);
    resultsFolds.RMSE = Inf(1,k);
    resultsFolds.RRMSE = Inf(1,k);
    resultsFolds.R2 = -Inf(1,k);
else
    numModels = trainParamsEnsemble.numTrees;
    resultsFolds.MAE = Inf(numModels,k);
    resultsFolds.MSE = Inf(numModels,k);
    resultsFolds.RMSE = Inf(numModels,k);
    resultsFolds.RRMSE = Inf(numModels,k);
    resultsFolds.R2 = -Inf(numModels,k);
end

% perform training and testing k times
for i = 1 : k
    if verbose
        disp(['Fold #' num2str(i)]);
    end
    Xtr = zeros(n-sizes(k-i+1), d);
    Ytr = zeros(n-sizes(k-i+1), 1);
    currsize = 0;
    for j = 1 : k
        if k-i+1 ~= j
            idxtrain = ind(offsets(j):offsets(j)+sizes(j)-1);
            Xtr(currsize+1:currsize+1+sizes(j)-1, :) = X(idxtrain, :);
            Ytr(currsize+1:currsize+1+sizes(j)-1, 1) = Y(idxtrain, 1);
            currsize = currsize + sizes(j);
        end
    end
    idxtst = ind(offsets(k-i+1):offsets(k-i+1)+sizes(k-i+1)-1);
    Xtst = X(idxtst, :);
    Ytst = Y(idxtst, 1);
    model = m5pbuild(Xtr, Ytr, trainParams, isBinCat, trainParamsEnsemble, verbose);
    res = m5ptest(model, Xtst, Ytst);
    resultsFolds.MAE(:,i) = res.MAE;
    resultsFolds.MSE(:,i) = res.MSE;
    resultsFolds.RMSE(:,i) = res.RMSE;
    resultsFolds.RRMSE(:,i) = res.RRMSE;
    resultsFolds.R2(:,i) = res.R2;
    if isempty(trainParamsEnsemble)
        [resultsFolds.nRules(:,i), vars] = countrulesandlistvars(model);
        resultsFolds.nVars(:,i) = length(vars);
    end
end

resultsTotal.MAE = mean(resultsFolds.MAE,2);
resultsTotal.MSE = mean(resultsFolds.MSE,2);
resultsTotal.RMSE = mean(resultsFolds.RMSE,2);
resultsTotal.RRMSE = mean(resultsFolds.RRMSE,2);
resultsTotal.R2 = mean(resultsFolds.R2,2);
if isempty(trainParamsEnsemble)
    resultsTotal.nRules = mean(resultsFolds.nRules,2);
    resultsTotal.nVars = mean(resultsFolds.nVars,2);
end
return
