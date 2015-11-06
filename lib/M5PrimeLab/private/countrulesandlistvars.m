function [nRules, vars] = countrulesandlistvars(model)
% Called from printinfo (m5pbuild) and m5pcv.
% Counts all rules (equal to the number of leaf nodes) in the tree and
% lists all original variables (not synthetic ones that are automatically
% made).

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

[nRules, usedVars] = countRulesAndListVarsDo(model.tree);
vars = [];
for v = 1:length(model.binCat.binCat)
    for u = usedVars(:).'
        if any(model.binCat.varMap{v} == u)
            vars = union(vars, v);
            break;
        end
    end
end
return

function [nRules, vars] = countRulesAndListVarsDo(node)
% Counts all rules (equal to the number of leaf nodes) in the tree and
% lists all (synthetic) variables.
if node.interior
    [nRules, vars] = countRulesAndListVarsDo(node.left);
    [nR, uV] = countRulesAndListVarsDo(node.right);
    nRules = nRules + nR;
    vars = union(vars, uV);
    vars = union(vars, node.splitAttr);
    if isfield(node, 'modelAttrIdx')
        vars = union(vars, node.modelAttrIdx);
    end
else
    nRules = 1;
    if isfield(node, 'modelAttrIdx')
        vars = node.modelAttrIdx;
    else
        vars = [];
    end
end
return
