function varImportance = aresanova(model, Xtr, Ytr, weights)
% aresanova
% Performs ANOVA decomposition and variable importance assessment of the
% given ARES model and reports the results. For details, see user's manual
% as well as Sections 3.5 and 4.3 in (Friedman, 1991a) and Sections 2.4,
% 4.1, and 4.4 in (Friedman, 1991b).
% The function works with single-response models only.
%
% Call:
%   varImportance = aresanova(model, Xtr, Ytr, weights)
%
% All the input arguments, except the last one, are required.
%
% Input:
%   model         : ARES model.
%   Xtr, Ytr      : Training data observations.
%   weights       : A vector of weights for observations. The same weights
%                   that were used when the model was built.
%
% Output:
%   varImportance : Relative variable importances. Scaled so that the
%                   relative importance of the most important variable has
%                   a value of 100.

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
if isempty(Xtr) || isempty(Ytr)
    error('Data is empty.');
end
n = size(Xtr,1);
if size(Ytr,1) ~= n
    error('The number of rows in Xtr and Ytr should be equal.');
end
if length(model) > 1
    error('This function works with single-response models only.');
else
    if iscell(model)
        model = model{1};
    end
end
if size(Ytr,2) ~= 1
    error('Ytr should have one column.');
end
if (nargin < 4)
    weights = [];
else
    if (~isempty(weights)) && ...
       ((size(weights,1) ~= n) || (size(weights,2) ~= 1))
        error('weights vector is of wrong size.');
    end
end

YtrVar = var2(Ytr, weights);

if model.trainParams.cubic
    fprintf('Type: piecewise-cubic\n');
else
    fprintf('Type: piecewise-linear\n');
end
fprintf('GCV: %g\n', model.GCV);
fprintf('R2GCV: %g\n', 1 - model.GCV / YtrVar);
fprintf('Total number of basis functions: %d\n', length(model.coefs));
fprintf('Total effective number of parameters: %g\n', ...
        length(model.coefs) + model.trainParams.c * length(model.knotdims) / 2);
fprintf('ANOVA decomposition:\n');
fprintf('Function\t\tSTD\t\t\t\tGCV\t\t\tR2GCV\t\t#basis\t#params\t\tvariable(s)\n');
counter = 0;
for i = 1 : model.trainParams.maxInteractions
    combs = nchoosek(1:length(model.minX),i);
    for j = 1 : size(combs,1)
        [modelReduced, usedBasis] = aresanovareduce(model, combs(j,:), true);
        if ~isempty(usedBasis)
            counter = counter + 1;
            fprintf('%d\t', counter);
            % standard deviation of the ANOVA function
            fprintf('%15f', sqrt(var2(arespredict(modelReduced, Xtr), weights)));
            % GCV when the basis functions of the ANOVA function are deleted
            modelReduced = aresdel(model, usedBasis, Xtr, Ytr, weights);
            fprintf('%16f', modelReduced.GCV);
            fprintf('%14.4f', 1 -  modelReduced.GCV / YtrVar);
            % the number of basis functions for that ANOVA function
            fprintf('%13d', length(usedBasis));
            % effective parameters
            fprintf('%9.2f\t\t', length(usedBasis) + model.trainParams.c * length(usedBasis) / 2);
            % used variables
            fprintf('%d ', combs(j,:));
            fprintf('\n');
        end
    end
end
fprintf('Relative variable importance:\n');
fprintf('Variable\tImportance\n');
nVars = length(model.minX);
nBasis = length(model.knotdims);
varImportance = zeros(nVars,1);
for v = 1 : nVars
    funcsToDel = [];
    for i = 1 : nBasis
        dims = model.knotdims{i};
        if any(dims == v)
            funcsToDel = [funcsToDel i];
        end
    end
    if isempty(funcsToDel)
        varImportance(v) = 0;
    else
        modelReduced = aresdel(model, funcsToDel, Xtr, Ytr, weights);
        varImportance(v) = sqrt(modelReduced.GCV) - sqrt(model.GCV);
    end
end
maxImp = max(varImportance);
varImportance = varImportance ./ maxImp .* 100;
for v = 1 : nVars
    fprintf('%d\t\t\t%10.3f\n', v, varImportance(v));
end
return

function res = var2(values, weights)
if isempty(weights)
    res = var(values, 1);
else
    valuesMean = sum(values(:,1) .* weights) / sum(weights);
    res = sum(((values(:,1) - valuesMean) .^ 2) .* weights) / sum(weights);
end
return
