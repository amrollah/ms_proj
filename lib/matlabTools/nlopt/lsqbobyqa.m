function [x, fopt, retcode] = lsqbobyqa(f, y, x0, xmin, xmax, xstep, opt, varargin)
%[x, fopt, retcode] = lsqbobyqa(f, y, x0, xmin, xmax, xstep, opt, varargin) 
%  f can be 
%    a string, then the vector to optimize must be named x
%    or a function handle f(x,varargin)
%  opt can be numeric, then it defines the max # function evaluations
%    (default = 1000)
%  varargin contains the arguments for the function to fit
%
% Example: x=lsqbobyqa('x(1)+x(2)*xx',yy,[1 1],[0 0],[2 2],[],[],xx)
%   does a linear regression for the data (xx,yy) with coefficients
%   between 0 and 2, starting with 1


if nargin<6, xstep=[]; end
if nargin<7, opt=[]; end
global lsqbobyqa_data
lsqbobyqa_data=[];
lsqbobyqa_data.f = f;
lsqbobyqa_data.y = y;
if ischar(f)
  for i = 8:nargin
    if isempty(inputname(i)), error('arguments for least squares fit must be identified by names'); end
    lsqbobyqa_data.d.(inputname(i)) = varargin{i-7};
  end
else
  lsqbobyqa_data.varargin = varargin;
end
if isempty(opt), opt=1000; end
if isnumeric(opt); maxeval=opt; opt=[]; opt.maxeval=maxeval; end
sx0 = size(x0);
x0 = x0(:)'; xmin = xmin(:)'; xmax = xmax(:)'; xstep = xstep(:)';
if isempty(xstep), xstep = (xmax-xmin)/20; end
opt.algorithm = NLOPT_LN_BOBYQA;
opt.lower_bounds = mop(0*x0,'+',xmin);
opt.upper_bounds = mop(0*x0,'+',xmax);
opt.initial_step = mop(0*x0,'+',xstep);
if ~isfield(opt,'xtol_rel'), opt.xtol_rel = 1e-8; end
opt.min_objective=@lsqbobyqa_ffit;
[x, fopt, retcode] = nloptmin(opt, x0);
x = reshape(x,sx0);

function y = lsqbobyqa_ffit(x)
global lsqbobyqa_data
if isfield(lsqbobyqa_data,'varargin')
  fev = feval(lsqbobyqa_data.f,x,lsqbobyqa_data.varargin{:});
else
  for lsqbobyqa_fn=fieldnames(lsqbobyqa_data.d)'
    eval([lsqbobyqa_fn{1} '=lsqbobyqa_data.d.' lsqbobyqa_fn{1} ';']);
  end
  fev = eval(lsqbobyqa_data.f);
end
y = sum((fev(:)-lsqbobyqa_data.y(:)).^2);
