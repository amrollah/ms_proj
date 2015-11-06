function Yq = predictsingle(node, Xq, modelTree)
% Called from m5ppredict and m5pbuild.
% Predicts response value for a single query point Xq.
% Assumes that Xq is already pre-processed, i.e., categorical variables
% replaced by binary in the same way as m5pbuild does it.

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

if node.interior
    % NaN is treated as the average of the corresponding
    % variable of the training observations reaching the node
    if (isnan(Xq(node.splitAttr)) && node.nanLeft) || ...
       (Xq(node.splitAttr) <= node.splitLocation)
        Yq = predictsingle(node.left, Xq, modelTree);
    else
        Yq = predictsingle(node.right, Xq, modelTree);
    end
else
    if modelTree
        if ~isempty(node.modelAttrIdx)
            % Replace NaNs with the average values of the corresponding
            % variables of the training observations reaching the node
            A = Xq(node.modelAttrIdx);
            where = isnan(A);
            A(where) = node.modelAttrAvg(where);
            % Calculate prediction
            Yq = [1 A] * node.modelCoefs;
        else
            Yq = node.modelCoefs;
        end
    else
        Yq = node.value;
    end
end
return
