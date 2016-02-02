function [y,g] = f1(x)
y = sum((x-0.3).^2,1);
if nargout>1, g = 2*(x-0.3); end;
