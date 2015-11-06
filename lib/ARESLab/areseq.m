function eq = areseq(model, precision, varNames, hideCubicSmoothing)
% areseq
% Prints ARES model.
% The function works with single-response models only.
%
% Call:
%   eq = areseq(model, precision, varNames)
%
% All the input arguments, except the first one, are optional.
% Piecewise-cubic spline equations are very complex as they involve
% smoothing. For easier interpretation use piecewise-linear models or set
% hideCubicSmoothing to true.
%
% Input:
%   model         : ARES model.
%   precision     : Number of digits in the model coefficients and knot
%                   sites. (default value = 15)
%   varNames      : A cell array of variable names to show instead of the
%                   generic ones. This is used for
%                   piecewise-linear models only.
%   hideCubicSmoothing : This is for piecewise-cubic models only. Whether
%                   to hide piecewise-cubic spline smoothing. It's easier
%                   to understand the equations if smoothing is hidden. The
%                   model will look like piecewise-linear but the
%                   coefficients will be for the actual piecewise-cubic
%                   model. (default value = true)
%
% Output:
%   eq            : A cell array of strings containing equations for
%                   individual basis functions and the main model.

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
if (nargin < 2) || isempty(precision)
    precision = 15;
end
if (nargin >= 3)
    if (~isempty(varNames)) && (length(varNames) ~= length(model.minX))
        error('Wrong number of names in varNames.');
    end
else
    varNames = [];
end
if (nargin < 4) || isempty(hideCubicSmoothing)
    hideCubicSmoothing = true;
end

p = ['%.' num2str(precision) 'g'];
eq = {};

% compose the individual basis functions
for i = 1 : length(model.knotdims)
    func = ['BF' num2str(i) ' ='];
    if model.parents(i) > 0
        func = [func ' BF' num2str(model.parents(i)) ' *'];
        start = length(model.knotdims{i});
    else
        start = 1;
    end
    for j = start : length(model.knotdims{i})
        
        if (~model.trainParams.cubic) || (hideCubicSmoothing)
            if model.knotdirs{i}(j) > 0 % here the hinge function looks like "_/"
                if model.knotsites{i}(j) > 0
                    m = '';
                else
                    m = '+';
                end
                if isempty(varNames)
                    func = [func ' max(0, x' num2str(model.knotdims{i}(j),p) ' ' m num2str(-model.knotsites{i}(j),p) ')'];
                else
                    func = [func ' max(0, ' varNames{model.knotdims{i}(j)} ' ' m num2str(-model.knotsites{i}(j),p) ')'];
                end
            else % here the hinge function looks like "\_"
                if isempty(varNames)
                    func = [func ' max(0, ' num2str(model.knotsites{i}(j),p) ' -x' num2str(model.knotdims{i}(j),p) ')'];
                else
                    func = [func ' max(0, ' num2str(model.knotsites{i}(j),p) ' -' varNames{model.knotdims{i}(j)} ')'];
                end
            end
        else
            t = num2str(model.knotsites{i}(j),p);
            t1 = num2str(model.t1(i,model.knotdims{i}(j)),p);
            t2 = num2str(model.t2(i,model.knotdims{i}(j)),p);
            pp = ['p' num2str(i) '_' num2str(j)];
            rr = ['r' num2str(i) '_' num2str(j)];
            d = num2str(model.knotdims{i}(j),p);
            f = ['f' num2str(i) '_' num2str(j)];
            if model.knotdirs{i}(j) > 0 % here the hinge function looks like "_/"
                iff = ['if (x' d ' <= ' t1 ') then ' f ' = 0'];
                disp(iff);
                eq{end+1,1} = iff;
                
                iff = ['if (' t1 ' < x' d ' < ' t2 ') then begin'];
                disp(iff);
                eq{end+1,1} = iff;
                pPoz = ['  ' pp ' = (2*(' t2 ') + (' t1 ') - 3*(' t ')) / ((' t2 ') - (' t1 '))^2'];
                disp(pPoz);
                eq{end+1,1} = pPoz;
                rPoz = ['  ' rr ' = (2*(' t ') - (' t2 ') - (' t1 ')) / ((' t2 ') - (' t1 '))^3'];
                disp(rPoz);
                eq{end+1,1} = rPoz;
                iff = ['  ' f ' = ' pp '*(x' d '-(' t1 '))^2 + ' rr '*(x' d '-(' t1 '))^3'];
                disp(iff);
                eq{end+1,1} = iff;
                iff = 'end';
                disp(iff);
                eq{end+1,1} = iff;
                
                iff = ['if (x' d ' >= (' t2 ')) then ' f ' = x' d ' - (' t ')'];
                disp(iff);
                eq{end+1,1} = iff;
                
                func = [func ' '  f];
            else % here the hinge function looks like "\_"
                iff = ['if (x' d ' <= ' t1 ') then ' f ' = -(x' d ' - (' t '))'];
                disp(iff);
                eq{end+1,1} = iff;
                
                iff = ['if (' t1 ' < x' d ' < ' t2 ') then begin'];
                disp(iff);
                eq{end+1,1} = iff;
                pNeg = ['  ' pp ' = (3*(' t ') - 2*(' t1 ') - (' t2 ')) / ((' t1 ') - (' t2 '))^2'];
                disp(pNeg);
                eq{end+1,1} = pNeg;
                rNeg = ['  ' rr ' = ((' t1 ') + (' t2 ') - 2*(' t ')) / ((' t1 ') - (' t2 '))^3'];
                disp(rNeg);
                eq{end+1,1} = rNeg;
                iff = ['  ' f ' = ' pp '*(x' d '-(' t2 '))^2 + ' rr '*(x' d '-(' t2 '))^3'];
                disp(iff);
                eq{end+1,1} = iff;
                iff = 'end';
                disp(iff);
                eq{end+1,1} = iff;
                
                iff = ['if (x' d ' >= (' t2 ')) then ' f ' = 0'];
                disp(iff);
                eq{end+1,1} = iff;
                
                func = [func ' '  f];
            end
        end
        
        if j < length(model.knotdims{i})
            func = [func ' *'];
        end
    end
    disp(func);
    eq{end+1,1} = func;
end

% compose the summation
func = ['y = ' num2str(model.coefs(1),p)];
for i = 1 : length(model.knotdims)
    if model.coefs(i+1) >= 0
        func = [func ' +'];
    else
        func = [func ' '];
    end
    func = [func num2str(model.coefs(i+1),p) '*BF' num2str(i)];
end
disp(func);
eq{end+1,1} = func;

if model.trainParams.cubic
    if hideCubicSmoothing
        fprintf('\nWARNING: Piecewise-cubic spline smoothing hidden.\n');
    else
        fprintf('\nINFO: Piecewise-cubic spline equations are very complex as they involve smoothing. For easier interpretation use piecewise-linear models or set hideCubicSmoothing to true.\n');
    end
end

return
