function trainParams = m5pparams(modelTree, minNumCases, prune, ...
    smoothingK, splitThreshold, aggressivePruning)
% m5pparams
% Creates configuration for building M5' trees. The structure is for
% further use with m5pbuild and m5pcv functions.
%
% Call:
%   trainParams = m5pparams(modelTree, minNumCases, prune, ...
%                           smoothingK, splitThreshold, aggressivePruning)
%
% All the input arguments of this function are optional. Empty values are
% also accepted (the corresponding default values will be used).
% It is quite possible that the default values for minNumCases and
% smoothingK will be far from optimal for your data.
% For a typical configuration for trees in Random Forests, set modelTree =
% false, minNumCases = 5, prune = false, smoothingK = 0 with all the other
% parameters left to their defaults.
%
% Input:
%   modelTree     : Whether to build a model tree (true) or a regression
%                   tree (false) (default value = true). For model trees,
%                   each leaf node will contain a linear regression model.
%                   If pruning is enabled, all the models will be reduced
%                   in size using sequential backward selection algorithm
%                   by greedily dropping terms.
%   minNumCases   : The minimum number of training observations one node
%                   may represent. Values lower than 2 are not allowed
%                   (default value = 4). If built trees are too large or
%                   overfit the data, consider increasing the number - this
%                   will result in smaller trees that are less sensitive
%                   to noise (but can also be underfitted).
%   prune         : Whether to prune the tree (default value = true).
%                   First, all models are reduced in size by greedily
%                   dropping terms and then the tree itself is pruned by
%                   eliminating leaves and subtrees if doing so improves
%                   the estimated error.
%   smoothingK    : Smoothing parameter. Larger values mean more smoothing.
%                   Usually not recommended for regression trees but can be
%                   useful for model trees. Set to 0 to disable. Default
%                   value = 15 (Quinlan, 1992; Wang & Witten, 1997).
%                   Smoothing tries to compensate for sharp discontinuities
%                   occurring between adjacent nodes of the tree. In case
%                   studies by Quinlan, 1992, as well as Wang & Witten,
%                   1997, this almost always had a positive effect on
%                   results. Smoothing is performed after building and
%                   pruning therefore this parameter does not influence
%                   those processes. Unfortunately, smoothed trees are
%                   harder to interpret.
%   splitThreshold : A node is not split if the standard deviation of the
%                   values of output variable at the node is less than
%                   splitThreshold of the standard deviation of response
%                   variable for the entire training data (default value
%                   = 0.05 (i.e., 5%) (Wang & Witten, 1997)). The results
%                   are usually not very sensitive to the exact choice of
%                   the threshold (Wang & Witten, 1997).
%   aggressivePruning : By default pruning is done as proposed by Quinlan,
%                   1992, and Wang & Witten, 1997, but you can also employ
%                   more aggressive pruning, the one that is implemented in
%                   Weka's version of M5' (Hall et al., 2009). Simply put,
%                   in the aggressive version, while estimating error of a
%                   subtree, one penalizes not only the number of
%                   parameters of regression models at its leaves but also
%                   its total number of splits. Aggressive pruning will
%                   produce significantly smaller trees that are easier to
%                   interpret but it can also cause underfitting. (default
%                   value = false)
%
% Output:
%   trainParams   : A structure of parameters for further use with m5pbuild
%                   and m5pcv functions containing the provided values (or
%                   default ones, if not provided).

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

if (nargin < 1) || isempty(modelTree)
    trainParams.modelTree = true;
else
    trainParams.modelTree = modelTree;
end

if (nargin < 2) || isempty(minNumCases)
    trainParams.minNumCases = 4;
else
    trainParams.minNumCases = max(2, minNumCases);
end

if (nargin < 3) || isempty(prune)
    trainParams.prune = true;
else
    trainParams.prune = prune;
end

if (nargin < 4) || isempty(smoothingK)
    trainParams.smoothingK = 15;
else
    trainParams.smoothingK = max(0, smoothingK);
end

if (nargin < 5) || isempty(splitThreshold)
    trainParams.splitThreshold = 0.05;
else
    trainParams.splitThreshold = max(0, splitThreshold);
end

if (nargin < 6) || isempty(aggressivePruning)
    trainParams.aggressivePruning = false;
else
    trainParams.aggressivePruning = aggressivePruning;
end

return
