function x = mycmaes(f,x0,lb,ub,sigma,print_iter,varargin)
x = x0;
opts=cmaes('defaults');
if nargin<3, lb = []; end
if nargin<4, ub = []; end
if nargin<5, sigma = []; end
if size(x0,2)>1, 
  x0 = x0'; lb = lb'; ub = ub';
end
if nargin<6, print_iter = 0; end
if ~print_iter
  opts.DispModulo = 0; opts.DispFinal = 'off';
end
opts.SaveVariables = 'off'; opts.LogModulo = 0; opts.LogTime = 0;
for i = 1:2:length(varargin)-1
  opts.(varargin{i}) = varargin{i+1};
end
x(:) = (lb+ub)/2;
if isempty(lb) || isempty(ub), j=true(size(x0));
else
  j = (ub>=lb+1e-8);
  x0 = x0(j); lb = lb(j); ub = ub(j);
end
if ~isempty(lb), opts.LBounds = lb; end
if ~isempty(ub), opts.UBounds = ub; end
if isempty(sigma)
  if ~isempty(lb) && ~isempty(ub), sigma=(ub-lb)/5;
  else sigma=ones(size(x0)); end
end
x(j) = cmaes(@fobj,x0,sigma,opts);
function y = fobj(x1)
  x(j) = x1;
  y = feval(f,x);
end
end
