function x = softunique(x,tol,varargin)
if size(x,1)==1, x=x'; end
x = unique(x,varargin{:});
j = max(abs(diff(x)),[],2)<tol;
j = [j; false];
x(j,:) = [];
if size(x,1)==1, x=x'; end
