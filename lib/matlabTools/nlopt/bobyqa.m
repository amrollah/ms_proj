function [x_sol, fopt, retcode] = bobyqa(f, x0, xmin, xmax, xstep, opt)
%[x, fopt, retcode] = bobyqa(f, x0, xmin, xmax, xstep)
%[x, fopt, retcode] = bobyqa(f, x0)
%  uses xmin=zeros(d,1), xmax=ones(d,1), xstep=0,02*(xmax-xmin)
%[x, fopt, retcode] = bobyqa(f, x0, xmin, xmax, xstep, opt)
%  opt can be numeric, then it defines the max # function evaluations

x_sol = x0;
x0 = x0(:)';
d = length(x0);
if nargin<3, xmin=[]; end
if nargin<4, xmax=[]; end
if nargin<5, xstep=[]; end
if nargin<6, opt=[]; end
if isempty(xmin), xmin=zeros(d,1); end
if isempty(xmax), xmax=ones(d,1); end
if isempty(xstep), xstep=0.02*(xmax-xmin); end
if isempty(opt), opt=1000; end
if isnumeric(opt); maxeval=opt; opt=[]; opt.maxeval=maxeval; end
xmin = xmin(:)'; xmax = xmax(:)'; xstep = xstep(:)';
j2opt = xmin<xmax;
x_sol(~j2opt) = (xmin(~j2opt)+xmax(~j2opt))/2;
x0 = x0(j2opt); xmin = xmin(j2opt); 
xmax = xmax(j2opt); xstep = xstep(j2opt);
x0 = min(xmax,max(xmin,x0));
opt.algorithm = 34; %NLOPT_LN_BOBYQA;
opt.lower_bounds = 0*x0+xmin;
opt.upper_bounds = 0*x0+xmax;
opt.initial_step = 0*x0+xstep;
if ~isfield(opt,'xtol_rel'), opt.xtol_rel = 1e-8; end
opt.min_objective=@bobyqa_fun;
[x_sol(j2opt), fopt, retcode] = nloptmin(opt, x0);

  function y = bobyqa_fun(x)
    x1 = x_sol;
    x1(j2opt) = x;
    y = feval(f,x1);
  end

end
