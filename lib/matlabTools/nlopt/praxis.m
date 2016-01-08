function [x_sol, fopt, retcode] = praxis(f, x0, xstep, opt)
%[x, fopt, retcode] = praxis(f, x0, xstep [,opt])
%  opt can be numeric, then it defines the max # function evaluations

x0 = x0(:)';
if nargin<4, opt=[]; end
if isempty(xstep), xstep=0.02*(xmax-xmin); end
if isempty(opt), opt=1000; end
if isnumeric(opt); maxeval=opt; opt=[]; opt.maxeval=maxeval; end
xstep = xstep(:)';
opt.algorithm = NLOPT_LN_PRAXIS;
opt.initial_step = 0*x0+xstep;
if ~isfield(opt,'xtol_rel'), opt.xtol_rel = 1e-8; end
opt.min_objective=f;
[x_sol, fopt, retcode] = nloptmin(opt, x0);
end
