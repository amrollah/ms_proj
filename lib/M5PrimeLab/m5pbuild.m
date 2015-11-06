function [model, time, ensembleResults] = m5pbuild(Xtr, Ytr, trainParams, ...
    isBinCat, trainParamsEnsemble, verbose)
% m5pbuild
% Builds M5' regression tree, model tree, or ensemble of trees.
%
% Call:
%   [model, time, ensembleResults] =
%       m5pbuild(Xtr, Ytr, trainParams, isBinCat, trainParamsEnsemble, verbose)
%
% All the input arguments, except the first two, are optional. Empty values
% are also accepted (the corresponding default values will be used).
%
% Input:
%   Xtr, Ytr      : Xtr is a matrix with rows corresponding to
%                   observations, and columns corresponding to input
%                   variables. Ytr is a column vector of response values.
%                   Input variables can be continuous, binary, as well as
%                   categorical (indicate using isBinCat). All values must
%                   be numeric. Missing values in Xtr must be indicated as
%                   NaN.
%   trainParams   : A structure of training parameters for the algorithm.
%                   If not provided, default values will be used (see
%                   function m5pparams for details).
%   isBinCat      : A vector of flags indicating type of each input
%                   variable - either continuous (false) or categorical
%                   with any number of categories, including binary (true).
%                   The vector should be of the same length as the number
%                   of columns in Xtr. m5pbuild then detects all the
%                   actually possible values for categorical variables from
%                   the training data. Any new values detected later, e.g.,
%                   in the test data, will be treated as NaN. By default,
%                   the vector is created with all values equal to false,
%                   meaning that all the variables are treated as
%                   continuous.
%   trainParamsEnsemble : A structure of parameters for building ensemble
%                   of trees. If not provided, a single tree is built. See
%                   function m5pparamsensemble for details. This can also
%                   be useful for variable importance assessment. See
%                   user's manual for examples of usage.
%                   Note that the ensemble building algorithm employs
%                   random number generator for which you can set seed
%                   before calling m5pbuild.
%   verbose       : Whether to output additional information to console.
%                   (default value = true)
%
% Output:
%   model         : A single M5' tree or a cell array of M5' trees if an
%                   ensemble is built. A structure defining one tree has
%                   the following fields:
%     binCat      : Information regarding original (continuous / binary /
%                   categorical) variables, transformed (synthetic binary)
%                   variables, possible values for categorical variables,
%                   and lowest values for all the variables.
%     trainParams : A structure of training parameters for the algorithm
%                   (updated if any values are chosen automatically).
%     tree        : A structure defining the built M5' tree.
%   time          : Algorithm execution time (in seconds).
%   ensembleResults : A structure of results if ensemble of tree is built.
%                   The structures has the following fields:
%     OOBError    : Out-of-bag estimate of prediction Mean Squared Error of
%                   the ensemble after each new tree is built. The number
%                   of rows is equal to the number of trees built. OOBError
%                   is available only if getOOBError in trainParamsEnsemble
%                   is set to true.
%     OOBNum      : Number of times observations were out-of-bag (and thus
%                   used in computing OOBError). The length of OOBNum is
%                   equal to the number of rows in Xtr and Ytr. OOBNum is
%                   available only if getOOBError in trainParamsEnsemble is
%                   set to true.
%     varImportance : Variable importance assessment. Calculated when
%                   out-of-bag data of a variable is permuted. A matrix
%                   with four rows and as many columns as there are columns
%                   in Xtr. First row is the average increase of out-of-bag
%                   Mean Absolute Error (MAE), second row is standard
%                   deviation of the average increase of MAE, third row is
%                   the average increase of out-of-bag Mean Squared Error
%                   (MSE), fourth row isstandard deviation of the average
%                   increase of MSE. The final variable importance estimate
%                   is usually calculated by dividing each MAE or MSE by
%                   the corresponding standard deviation. Bigger values
%                   then indicate bigger importance of the corresponding
%                   variable. See user's manual for example of usage.
%                   varImportance is available only if getVarImportance in
%                   trainParamsEnsemble is > 0.

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

% Citing the M5PrimeLab toolbox:
% Jekabsons G., M5PrimeLab: M5' regression tree, model tree, and tree ensemble
% toolbox for Matlab/Octave, 2015, available at http://www.cs.rtu.lv/jekabsons/

% Last update: October 6, 2015

if nargin < 2
    error('Not enough input arguments.');
end

if isempty(Xtr) || isempty(Ytr)
    error('Training data is empty.');
end
[n, mOriginal] = size(Xtr); % number of observations and number of input variables
if size(Ytr,1) ~= n
    error('The number of rows in Xtr and Ytr should be equal.');
end
if size(Ytr,2) ~= 1
    error('Ytr should have one column.');
end

if (nargin < 3) || isempty(trainParams)
    trainParams = m5pparams();
else
    trainParams.minNumCases = max(2, trainParams.minNumCases);
end
if (nargin < 4) || isempty(isBinCat)
    isBinCat = false(1,mOriginal);
else
    isBinCat = isBinCat(:).'; % force row vector
    if length(isBinCat) ~= mOriginal
        error('The number of elements in isBinCat should be equal to the number of columns in Xtr.');
    end
end
if (nargin < 5)
    trainParamsEnsemble = [];
end
if (nargin < 6) || isempty(verbose)
    verbose = true;
end

binCat = isBinCat .* 2;
% Fill binCat with numbers of unique values
for i = 1 : mOriginal
    if binCat(i) >= 2
        u = unique(Xtr(:,i));
        u = u(~isnan(u)); % no NaNs
        len = length(u);
        binCat(i) = len;
        if len >= 50
            disp(['Warning: Categorical variable #' num2str(i) ' has ' num2str(len) ' unique values.']);
        end
    end
end
% Transform categorical variables into a number of synthetic binary variables
% binCat < 2 for continuous variables
% binCat = 2 for original binary variables
% binCat > 2 for binary variables created from original categorical variables
binCatVals = {};
if any(binCat >= 2)
    binCatNew = [];
    binCatCounter = 0;
    Xnew = [];
    model.binCat.varMap = {};
    for i = 1 : mOriginal
        if binCat(i) >= 2
            XX = Xtr(:,i);
            u = unique(XX(~isnan(XX))); % no NaNs, sorted
            if size(u,1) > 2
                model.binCat.varMap = [model.binCat.varMap (size(binCatNew,2)+1) : (size(binCatNew,2)+size(u,1)-1)];
                avg = zeros(size(u,1),1);
                for j = 1 : size(u,1)
                    avg(j) = mean(Ytr(Xtr(:,i) == u(j)));
                end
                [dummy, ind] = sort(avg);
                u = u(ind);
                Xb = zeros(n,size(u,1)-1);
                for j = 1 : n
                    if isnan(Xtr(j,i))
                        Xb(j,:) = NaN;
                    else
                        Xb(j, 1 : find(Xtr(j,i) == u) - 1) = 1;
                    end
                end
                Xnew = [Xnew Xb];
                binCatNew = [binCatNew repmat(size(u,1),1,size(u,1)-1)];
            else
                Xnew = [Xnew Xtr(:,i)];
                binCatNew = [binCatNew 2];
                model.binCat.varMap = [model.binCat.varMap size(binCatNew,2)];
            end
            binCat(i) = size(u,1);
            binCatCounter = binCatCounter + 1;
            binCatVals{binCatCounter} = u;
        else
            Xnew = [Xnew Xtr(:,i)];
            binCatNew = [binCatNew binCat(i)];
            model.binCat.varMap = [model.binCat.varMap size(binCatNew,2)];
        end
    end
    Xtr = Xnew;
    model.binCat.catVals = binCatVals;
else
    binCatNew = binCat;
    model.binCat.varMap = num2cell(1:mOriginal);
end

model.binCat.binCat = binCat;
model.binCat.binCatNew = binCatNew >= 2; % 0 for continuous; 1 for binary
if any(model.binCat.binCatNew)
    % this is used only for proper output of the tree when there are binary variables (might be synthetic)
    model.binCat.minVals = min(Xtr);
end

model.trainParams = trainParams;

if verbose
    if trainParams.modelTree, str = 'model'; else str = 'regression'; end
    if isempty(trainParamsEnsemble)
        disp(['Growing M5'' ' str ' tree...']);
    else
        disp(['Growing M5'' ' str ' tree ensemble...']);
    end
end
ws = warning('off');
ttt = tic;

% For original binary and continuous variables beta = 1
% For synthetic binary variables created from original categorical variables beta < 1
beta = exp(7 * (2 - max(2, binCatNew)) / n);

if isempty(trainParamsEnsemble)
    
    sd = std2(Ytr);
    numNotMissing = sum(~isnan(Xtr),1); % number of non-missing values for each variable
    % Growing the tree
    model.tree = splitNode(Xtr, Ytr, 1:n, sd, numNotMissing, binCatNew, trainParams, beta, [], [], []);
    if verbose && trainParams.prune
        disp('Pruning...');
    end
    % Pruning the tree and filling it with models
    model.tree = pruneNode(model.tree, Xtr, Ytr, trainParams);
    if trainParams.smoothingK > 0
        totalAttrs = model.binCat.varMap{end}(end);
        model.tree = smoothing(model.tree, [], trainParams.modelTree, trainParams.smoothingK, totalAttrs);
    end
    model.tree = cleanUp(model.tree, trainParams.modelTree);
    
    ensembleResults = [];
    
else
    
    if trainParamsEnsemble.numVarsTry < 1
        if trainParamsEnsemble.numVarsTry < 0
            trainParamsEnsemble.numVarsTry = mOriginal / 3;
        else
            trainParamsEnsemble.numVarsTry = mOriginal;
        end
    end
    trainParamsEnsemble.numVarsTry = min(mOriginal, max(1, floor(trainParamsEnsemble.numVarsTry)));
    
    if round(trainParamsEnsemble.inBagFraction * n) < 1
        error('trainParamsEnsemble.inBagFraction too small. In-bag set empty.');
    end
    if (~trainParamsEnsemble.withReplacement) && (round(trainParamsEnsemble.inBagFraction * n) >= n)
        error('trainParamsEnsemble.inBagFraction too big. Out-of-bag set empty.');
    end
    
    modelBase = model;
    models = cell(trainParamsEnsemble.numTrees, 1);
    
    if (~trainParamsEnsemble.getOOBError) && (trainParamsEnsemble.getVarImportance == 0)
        ensembleResults = [];
    else
        if trainParamsEnsemble.getOOBError
            OOBPred = zeros(n, 1);
            OOBNum = zeros(n, 1);
            ensembleResults.OOBError = NaN(trainParamsEnsemble.numTrees, 1);
        end
        if trainParamsEnsemble.getVarImportance > 0
            diffOOBMAE = NaN(trainParamsEnsemble.numTrees, mOriginal);
            diffOOBMSE = NaN(trainParamsEnsemble.numTrees, mOriginal);
            ensembleResults.varImportance = zeros(4, mOriginal); % increase in MAE, SD, increase in MSE, SD
        end
    end
    
    % for each tree
    for t = 1 : trainParamsEnsemble.numTrees
        if verbose && (trainParamsEnsemble.verboseNumIter > 0) && ...
                (mod(t, trainParamsEnsemble.verboseNumIter) == 0)
            fprintf('Growing tree #%d...\n', t);
        end
        % sampling
        if trainParamsEnsemble.withReplacement
            idx = randi(n, round(trainParamsEnsemble.inBagFraction * n), 1);
            X = Xtr(idx,:);
            Y = Ytr(idx,1);
        else
            perm = randperm(n);
            idx = perm(1:round(trainParamsEnsemble.inBagFraction * n));
            X = Xtr(idx,:);
            Y = Ytr(idx,1);
        end
        
        if t > 1
            model = modelBase;
        end
        
        sd = std2(Y);
        numNotMissing = sum(~isnan(X),1); % number of non-missing values for each variable
        % Growing the tree
        model.tree = splitNode(X, Y, 1:size(Y,1), sd, numNotMissing, binCatNew, trainParams, beta, ...
            trainParamsEnsemble.numVarsTry, mOriginal, model.binCat.varMap);
        % Pruning the tree and filling it with models
        model.tree = pruneNode(model.tree, X, Y, trainParams);
        if trainParams.smoothingK > 0
            totalAttrs = model.binCat.varMap{end}(end);
            model.tree = smoothing(model.tree, [], trainParams.modelTree, trainParams.smoothingK, totalAttrs);
        end
        model.tree = cleanUp(model.tree, trainParams.modelTree);
        
        % additional calculations, if asked
        if trainParamsEnsemble.getOOBError || (trainParamsEnsemble.getVarImportance > 0)
            idxoob = true(n,1);
            idxoob(idx) = false;
            idxoob = find(idxoob);
            if ~isempty(idxoob) % test for the unlikely case when out-of-bag set is empty
                Xoob = Xtr(idxoob,:);
                Yq = zeros(size(Xoob,1),1);
                for i = 1 : size(Xoob,1)
                    Yq(i) = predictsingle(model.tree, Xoob(i,:), trainParams.modelTree);
                end
                if trainParamsEnsemble.getOOBError
                    OOBNum(idxoob) = OOBNum(idxoob) + 1;
                    idxNonzero = OOBNum ~= 0;
                    OOBPred(idxoob) = OOBPred(idxoob) + Yq;
                    ensembleResults.OOBError(t,1) = mean(((OOBPred(idxNonzero) ./ OOBNum(idxNonzero)) - Ytr(idxNonzero)) .^ 2);
                    if verbose && (trainParamsEnsemble.verboseNumIter > 0) && ...
                            (mod(t, trainParamsEnsemble.verboseNumIter) == 0)
                        fprintf('Out-of-bag MSE with %d trees: %.5g\n', t, ensembleResults.OOBError(t,1));
                    end
                end
                
                if trainParamsEnsemble.getVarImportance > 0
                    Yqtdiff = Yq - Ytr(idxoob);
                    for v = 1 : mOriginal
                        for iPerm = 1:trainParamsEnsemble.getVarImportance
                            Xoobpert = Xoob;
                            idxoobpert = idxoob(randperm(size(idxoob,1)),1);
                            % Perturb OOB variables that correspond to the original vth variable
                            for vnew = model.binCat.varMap{v}
                                Xoobpert(:,vnew) = Xtr(idxoobpert,vnew);
                            end
                            Yqpert = zeros(size(Xoobpert,1),1);
                            for i = 1 : size(Xoobpert,1)
                                Yqpert(i) = predictsingle(model.tree, Xoobpert(i,:), trainParams.modelTree);
                            end
                            Yqptdiff = Yqpert - Ytr(idxoob);
                            if iPerm == 1
                                diffOOBMAE(t,v) = mean(abs(Yqptdiff)) - mean(abs(Yqtdiff));
                                diffOOBMSE(t,v) = mean(Yqptdiff .^ 2) - mean(Yqtdiff .^ 2);
                            else
                                diffOOBMAE(t,v) = diffOOBMAE(t,v) + mean(abs(Yqptdiff)) - mean(abs(Yqtdiff));
                                diffOOBMSE(t,v) = diffOOBMSE(t,v) + mean(Yqptdiff .^ 2) - mean(Yqtdiff .^ 2);
                            end
                        end
                        if trainParamsEnsemble.getVarImportance > 1
                            diffOOBMAE(t,v) = diffOOBMAE(t,v) / trainParamsEnsemble.getVarImportance;
                            diffOOBMSE(t,v) = diffOOBMSE(t,v) / trainParamsEnsemble.getVarImportance;
                        end
                    end
                end
                
            end
        end
        
        models{t} = model;
    end % end of loop through all trees
    model = models;
    if trainParamsEnsemble.getOOBError
        ensembleResults.OOBNum = OOBNum;
    end
    if trainParamsEnsemble.getVarImportance > 0
        ensembleResults.varImportance(1,:) = mean(diffOOBMAE, 1);
        ensembleResults.varImportance(2,:) = std(diffOOBMAE, 1, 1);
        ensembleResults.varImportance(3,:) = mean(diffOOBMSE, 1);
        ensembleResults.varImportance(4,:) = std(diffOOBMSE, 1, 1);
    end
    
end

time = toc(ttt);
if verbose
    if isempty(trainParamsEnsemble)
        printinfo(model);
    end
    fprintf('Execution time: %0.2f seconds\n', time);
end
warning(ws);
return

%==========================================================================

function [node, attrList] = splitNode(X, Y, caseInd, sd, numNotMissing, binCat, trainParams, beta, numVarsTry, mOriginal, varMap)
% Splits node into left node and right node
node.caseInd = caseInd;
YY = Y(caseInd);
if (length(caseInd) < trainParams.minNumCases * 2) || ...
   (std(YY) < trainParams.splitThreshold * sd)
    node.interior = false; % this node will be a leaf node
    attrList = [];
    return
end;
sdr = -Inf;
if isempty(numVarsTry) || (numVarsTry >= mOriginal)
    varsTry = 1:size(X, 2); % we will try all variables
else
    % we will try random subset of variables
    % if categorical variable is chosen, we will try all its synthetic binary variables
    randm = randperm(mOriginal);
    varsTry = [];
    for v = randm(1:numVarsTry)
        for vnew = varMap{v}
            varsTry = [varsTry vnew];
        end
    end
end
% let's find best variable and best split
for i = varsTry
    XX = X(caseInd,i);
    % NaNs will not be used for split point determination
    % and there is no need to sort because unique already sorts
    sorted = unique(XX(~isnan(XX)));
    if size(sorted,1) < 2
        continue;
    end
    splitCandidates = (sorted(1:end-1) + sorted(2:end)) ./ 2;
    % let's find best split
    for j = 1 : size(splitCandidates,1)
        sdrtmp = splitsdr(splitCandidates(j), XX, YY, ...
                          numNotMissing(i), trainParams.minNumCases, beta(i));
        if sdrtmp > sdr
            sdr = sdrtmp;
            splitPoint = splitCandidates(j);
            attr = i;
        end
    end
end
if sdr <= 0
    node.interior = false; % this node will be a leaf node
    attrList = [];
else
    [leftInd, rightInd] = leftright(splitPoint, X(caseInd,attr), YY, binCat(attr));
    leftInd = caseInd(leftInd);
    rightInd = caseInd(rightInd);
    node.interior = true; % this node will be an interior node
    node.splitAttr = attr;
    node.splitLocation = splitPoint;
    attrList = attr;
    [node.left, attrList2] = splitNode(X, Y, leftInd, sd, numNotMissing, binCat, trainParams, beta, numVarsTry, mOriginal, varMap);
    attrList = [attrList attrList2];
    [node.right, attrList2] = splitNode(X, Y, rightInd, sd, numNotMissing, binCat, trainParams, beta, numVarsTry, mOriginal, varMap);
    attrList = unique([attrList attrList2]);
    node.attrList = attrList;
end
return

function sdr = splitsdr(split, X, Y, numNotMissing, minNumCases, beta)
% Calculates SDR for the specific split point
leftInd = find(X <= split);
if (size(leftInd,1) < minNumCases)
    sdr = -Inf;
    return
end
rightInd = find(X > split);
if (size(rightInd,1) < minNumCases)
    sdr = -Inf;
    return
end
sizeBoth = size(leftInd,1) + size(rightInd,1); % size without missing values
if size(Y,1) == sizeBoth % if all values for the variable are known
    sdr = std2(Y) - (size(leftInd,1) * std2(Y(leftInd)) + size(rightInd,1) * std2(Y(rightInd))) / sizeBoth;
else
    %{
    allInd = [leftInd; rightInd];
    sdr = std2(Y(allInd)) - ...
          (size(leftInd,1) * std2(Y(leftInd)) + size(rightInd,1) * std2(Y(rightInd))) / sizeBoth;
    %}
    % calculation without making a new vector
    sdr = std3(Y(leftInd), Y(rightInd)) - ...
          (size(leftInd,1) * std2(Y(leftInd)) + size(rightInd,1) * std2(Y(rightInd))) / sizeBoth;
end
sdr = sdr * numNotMissing / sizeBoth * beta;
return

function stdev = std2(Y)
% Calculates standard deviation
% Does the same as Matlab's std function but without all the overhead
nn = size(Y,1);
stdev = sqrt(sum((Y - (sum(Y) / nn)) .^ 2) / (nn - 1));
return

function stdev = std3(Y1, Y2)
% Calculates standard deviation for two vectors
nn = size(Y1,1) + size(Y2,1);
avgY = (sum(Y1) + sum(Y2)) / nn;
sumSq = sum((Y1 - avgY) .^ 2);
sumSq = sumSq + sum((Y2 - avgY) .^ 2);
stdev = sqrt(sumSq / (nn - 1));
return

function [leftInd, rightInd] = leftright(split, X, Y, binCat)
% Splits all observations into left and right sets. Deals with NaNs separately.
leftInd = find(X <= split);
rightInd = find(X > split);
% Place observations with NaNs in left or right according to their Y values
isNaN = isnan(X);
if any(isNaN)
    if binCat < 2
        % For continuous variables
        [dummy, sorted] = sort(X(leftInd));
        sorted = leftInd(sorted);
        leftAvg = mean(Y(sorted(end - min([2 size(leftInd,1)-1]) : end)));
        [dummy, sorted] = sort(X(rightInd));
        sorted = rightInd(sorted);
        rightAvg = mean(Y(sorted(1 : min([3 size(rightInd,1)]))));
    else
        % For both original and synthetic binary variables
        leftAvg = mean(Y(leftInd));
        rightAvg = mean(Y(rightInd));
    end
    avgAvg = (leftAvg + rightAvg) / 2;
    smaller = Y(isNaN) <= avgAvg;
    fn = find(isNaN);
    if leftAvg <= rightAvg
        leftInd = [leftInd; fn(smaller)];
        rightInd = [rightInd; fn(~smaller)];
    else
        leftInd = [leftInd; fn(~smaller)];
        rightInd = [rightInd; fn(smaller)];
    end
end
return

function node = pruneNode(node, X, Y, trainParams)
% Prunes the tree and fills it with models (or average values).
% If tree pruning is disabled, only filling with models is done.
% For each model, subset selection is done (using backward selection).
if ~node.interior
    if ~trainParams.modelTree
        node.value = mean(Y(node.caseInd));
    else
        node.modelCoefs = mean(Y(node.caseInd));
        node.modelAttrIdx = [];
    end
    return;
end
attrInd = node.attrList;
node.left = pruneNode(node.left, X, Y, trainParams);
node.right = pruneNode(node.right, X, Y, trainParams);
if ~trainParams.modelTree
    node.value = mean(Y(node.caseInd));
    errNode = calcErrNodeWithAllKnown(node, X, Y, trainParams, true); % pretend known because regression tree doesn't care
else
    if isempty(attrInd) % no attributes. model will have only intercept
        node.modelCoefs = mean(Y(node.caseInd));
        node.modelAttrIdx = [];
        errNode = calcErrNodeWithAllKnown(node, X, Y, trainParams, true); % pretend known because no attributes are used
    else
        XX = X;
        isNaN = isnan(X(node.caseInd,attrInd));
        for i = 1 : length(attrInd)
            % Store average values of the variables (required when the tree is
            % used for prediction and NaN is encountered)
            % (node.modelAttrIdx provides index for variable for which
            % modelAttrAvg is the average value)
            node.modelAttrAvg(i) = mean(X(node.caseInd(~isNaN(:,i)),attrInd(i)));
            % Replace NaNs by the average values of the corresponding variables
            % of the training observations reaching the node
            XX(node.caseInd(isNaN(:,i)),attrInd(i)) = node.modelAttrAvg(i);
        end
        % Perform variable selection
        A = [ones(length(node.caseInd),1) XX(node.caseInd,attrInd)];
        node.modelCoefs = A \ Y(node.caseInd);
        node.modelAttrIdx = attrInd;
        errNode = calcErrNodeWithAllKnown(node, XX, Y, trainParams, true);
        if trainParams.prune
            attrIndBest = attrInd;
            coefsBest = node.modelCoefs;
            changed = false;
            for j = 1 : length(attrInd)
                attrIndOld = node.modelAttrIdx;
                for i = 1 : length(attrIndOld)
                    node.modelAttrIdx = attrIndOld;
                    node.modelAttrIdx(i) = [];
                    A = [ones(length(node.caseInd),1) XX(node.caseInd,node.modelAttrIdx)];
                    node.modelCoefs = A \ Y(node.caseInd);
                    errTry = calcErrNodeWithAllKnown(node, XX, Y, trainParams, true);
                    if errTry < errNode
                        attrIndBest = node.modelAttrIdx;
                        errNode = errTry;
                        coefsBest = node.modelCoefs;
                        changed = true;
                    end
                end
                node.modelAttrIdx = attrIndBest;
                node.modelCoefs = coefsBest;
                if ~changed
                    break;
                end
            end
            % Update node.modelAttrAvg if the used subset of variables has changed
            if length(node.modelAttrIdx) < length(attrInd)
                for i = 1 : length(node.modelAttrIdx)
                    node.modelAttrAvg(i) = node.modelAttrAvg(attrInd == node.modelAttrIdx(i));
                end
                node.modelAttrAvg = node.modelAttrAvg(1:length(node.modelAttrIdx));
            end
        end
    end
end
if trainParams.prune && ...
   ( ...
    ((~trainParams.aggressivePruning) && (calcErrSubtree(node, X, Y, trainParams) >= errNode)) || ...
    (trainParams.aggressivePruning && (calcErrSubtreeAggressive(node, X, Y, trainParams) >= errNode)) ...
   )
    % above we could also add "(sd * 0.000001 > errNode)"
    % this node will be a leaf node
    node.interior = false;
    node = rmfield(node, {'splitAttr', 'splitLocation', 'left', 'right', 'attrList'});
else
    % Store average value of the split variable (required when the tree
    % is used for prediction and NaN is encountered)
    notNaN = node.caseInd(~isnan(X(node.caseInd,node.splitAttr)));
    %node.splitAttrAvg = mean(X(notNaN,node.splitAttr)); % not really needed. we can just set nanLeft
    node.nanLeft = mean(X(notNaN,node.splitAttr)) <= node.splitLocation;
end
return

function err = calcErrSubtree(node, X, Y, trainParams)
% Calculates error of the subtree
if node.interior
    err = (length(node.left.caseInd) * calcErrSubtree(node.left, X, Y, trainParams) + ...
           length(node.right.caseInd) * calcErrSubtree(node.right, X, Y, trainParams)) / ...
           length(node.caseInd);
else
    err = calcErrNode(node, X, Y, trainParams);
end
return

function err = calcErrSubtreeAggressive(node, X, Y, trainParams)
% Calculates error of the subtree, applies penalty
[err, v] = calcErrSubtreeAggressiveDo(node, X, Y, trainParams);
nn = length(node.caseInd);
if (nn > v)
    err = err * (nn + v * 2) / (nn - v);
else
    err = err * 10;
end
return
function [err, v] = calcErrSubtreeAggressiveDo(node, X, Y, trainParams)
% Calculates error of the subtree
if node.interior
    [errLeft, vLeft] = calcErrSubtreeAggressiveDo(node.left, X, Y, trainParams);
    [errRight, vRight] = calcErrSubtreeAggressiveDo(node.right, X, Y, trainParams);
    err = (length(node.left.caseInd) * errLeft + length(node.right.caseInd) * errRight) / ...
           length(node.caseInd);
    v = vLeft + vRight + 1;
else
    err = calcErrNode(node, X, Y, trainParams);
    if trainParams.modelTree
        v = length(node.modelCoefs);
    else
        v = 1;
    end
end
return

function err = calcErrNode(node, X, Y, trainParams)
% Calculates error of the node. Handles missing values.
if trainParams.modelTree
    % Replace NaNs by the average values of the corresponding variables of
    % the training observations reaching the node
    isNaN = isnan(X(node.caseInd,node.modelAttrIdx));
    for i = 1 : length(node.modelAttrIdx)
        X(node.caseInd(isNaN(:,i)),node.modelAttrIdx(i)) = node.modelAttrAvg(i);
    end
end
err = calcErrNodeWithAllKnown(node, X, Y, trainParams, false);
return

function err = calcErrNodeWithAllKnown(node, X, Y, trainParams, forDroppingTerms)
% Calculates error of the node. Assumes all values are known.
if trainParams.modelTree
    val = [ones(length(node.caseInd),1) X(node.caseInd,node.modelAttrIdx)] * node.modelCoefs;
    deviation = mean(abs(val - Y(node.caseInd)));
    v = length(node.modelCoefs);
else
    deviation = mean(abs(node.value - Y(node.caseInd)));
    v = 1;
end
if ~trainParams.aggressivePruning
    nn = length(node.caseInd);
    err = (nn + v) / (nn - v) * deviation;
else
    if forDroppingTerms
        nn = length(node.caseInd);
        err = (nn + v * 2) / (nn - v) * deviation;
    else
        err = deviation;
    end
end
return

function node = cleanUp(node, modelTree)
% Removing the temporary fields
node.numCases = length(node.caseInd);
node = rmfield(node, 'caseInd');
if node.interior
    node = rmfield(node, 'attrList');
    if modelTree
        node = rmfield(node, 'modelCoefs');
        node = rmfield(node, 'modelAttrAvg');
        node = rmfield(node, 'modelAttrIdx');
    else
        node = rmfield(node, 'value');
    end
    node.left = cleanUp(node.left, modelTree);
    node.right = cleanUp(node.right, modelTree);
end
return

function node = smoothing(node, list, modelTree, smoothingK, totalAttrs)
% Performs smoothing by incorporating interior models into leaf models.
% Deals with modelAttrAvg, so that unknown values can be substituted with
% modelAttrAvg at leaves.
if node.interior
    if modelTree
        data.attrIdx = node.modelAttrIdx;
        data.coefs = node.modelCoefs;
        data.attrAvg = zeros(totalAttrs,1);
        data.attrAvg(node.modelAttrIdx) = node.modelAttrAvg;
    else
        data.value = node.value;
    end
    data.numCases = length(node.caseInd);
    list{end+1} = data; % making a list. will be used at leaf nodes
    node.left = smoothing(node.left, list, modelTree, smoothingK, totalAttrs);
    node.right = smoothing(node.right, list, modelTree, smoothingK, totalAttrs);
else
    if modelTree
        len = length(list);
        if len > 0
            attrIdx = node.modelAttrIdx;
            s_n = length(node.caseInd);
            coefs = zeros(totalAttrs+1,1);
            coefs([1 attrIdx+1]) = node.modelCoefs;
            attrAvg = zeros(totalAttrs,1);
            if ~isempty(attrIdx)
                attrAvg(attrIdx) = node.modelAttrAvg;
            end
            % pretend to go from the leaf node to the root node
            for i = len:-1:1
                % Update list of used variables
                attrIdx = union(attrIdx, list{i}.attrIdx);
                % Coefs at this node
                coefsHere = zeros(size(coefs));
                coefsHere([1 list{i}.attrIdx+1]) = list{i}.coefs;
                % Recalculate weighted averages for NaNs
                idx = true(size(coefs));
                idx(1) = false;
                idx((coefs == 0) & (coefsHere == 0)) = false;
                idxAttr = idx(2:end);
                attrAvg(idxAttr) = ...
                    attrAvg(idxAttr) .* s_n .* coefs(idx) ./ (s_n .* coefs(idx) + smoothingK .* coefsHere(idx)) + ...
                    list{i}.attrAvg(idxAttr) .* smoothingK .* coefsHere(idx) ./ (s_n .* coefs(idx) + smoothingK .* coefsHere(idx));
                % Recalculate smoothed coefs
                coefs = (s_n * coefs + smoothingK * coefsHere) / (s_n + smoothingK);
                s_n = list{i}.numCases; % s_n for next iteration
            end
            [attrIdx, idx] = sort(attrIdx); % sort attributes so that equations are easier to understand
            attrIdx = attrIdx(:)';
            node.modelCoefs = coefs([1 attrIdx(idx)+1]);
            node.modelAttrIdx = attrIdx;
            node.modelAttrAvg = attrAvg(attrIdx(idx))';
        end
    else
        len = length(list);
        if len > 0
            value = node.value;
            s_n = length(node.caseInd);
            % pretend to go from the leaf node to the root node
            for i = len:-1:1
                value = (s_n * value + smoothingK * list{i}.value) / (s_n + smoothingK); % calculate smoothed values
                s_n = list{i}.numCases; % s_n for next iteration
            end
            node.value = value;
        end
    end
end
return
