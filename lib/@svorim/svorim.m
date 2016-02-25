function m = svorim(x,y,kernel,c)
% m = p5model(x,y,kernel,c)
%   x: d*nx matrix
%   y: 1*nx vector
%   kernel: optional, either numeric (degree of polynomial) or struct 
%           defining n and b (polynomial) or struct defining sig
%           (Gaussian), default = 2
%   c: optional, default depends on the kernel
%   omitsmall: optional, default = 1
%   omitbig: optional, default = 1
%   options: optional options for ipopt

if nargin<3, kernel = []; end
if nargin<4, c=[]; end

if isempty(kernel), kernel = 2; end
if isempty(c), c=1000; end

[d,n] = size(x);
[y,jj] = sort(y);
x = x(:,jj);
  
if isnumeric(kernel)
  kerneln = kernel;
  kernel = [];
  kernel.n = kerneln;
end

if isfield(kernel,'sig')
  kerneltype = 0;
  kernelp1 = (d/kernel.sig)^2;
  kernelp2 = 0;
else
  kerneltype = 1;
  if ~isfield(kernel,'n'), kernel.n=1; end
  if ~isfield(kernel,'b'), kernel.b = 1; end
  if kernel.n<=1, kernel.b = 0; end
  kernelp1 = kernel.n;
  kernelp2 = kernel.b;
end

m.x = x;
m.kernel = kernel;
m.c = c;
a = svorim_train_mex(x,y,kerneltype,kernelp1,kernelp2,c);
m.a = sum(a.*(1-2*(repmat(y,size(a,1),1)<=repmat((1:size(a,1))',1,n))),1);
m = class(m,'svorim');
