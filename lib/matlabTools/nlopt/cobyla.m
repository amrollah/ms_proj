function [x, fopt, retcode] = cobyla(f, x0, xmin, xmax, xstep, opt, varargin)
%[x, fopt, retcode] = cobyla(f, x0, xmin, xmax, xstep)
%[x, fopt, retcode] = cobyla(f, x0)
%  uses xmin=zeros(d,1), xmax=ones(d,1), xstep=0,02*(xmax-xmin)
%[x, fopt, retcode] = cobyla(f, x0, xmin, xmax, xstep, opt) 
%  opt can be numeric, then it defines the max # function evaluations
%[x, fopt, retcode] = cobyla(f, x0, xmin, xmax, xstep, fc1) ==> rhs is 0
%[x, fopt, retcode] = cobyla(f, x0, xmin, xmax, xstep, fc1, rhs1, fc2, ...)
%[x, fopt, retcode] = cobyla(f, x0, xmin, xmax, xstep, opt, fc1, fc2, ...)
%etc.

constr = varargin;
[d,n] = size(x0);
if nargin<3, xmin=[]; end
if nargin<4, xmax=[]; end
if nargin<5, xstep=[]; end
if nargin<6, opt=[];
elseif ~isstruct(opt) && ~isempty(opt) && ~isnumeric(opt)
  constr = [{opt} constr];
  opt = [];
end
if isempty(xmin), xmin=zeros(d,1); end
if isempty(xmax), xmax=ones(d,1); end
if isempty(xstep), xstep=0.02*(xmax-xmin); end
if isempty(opt), opt=1000; end
if isnumeric(opt); maxeval=opt; opt=[]; opt.maxeval=maxeval; end
x0 = x0(:)'; xmin = xmin(:)'; xmax = xmax(:)'; xstep = xstep(:)';
opt.algorithm = 25; %NLOPT_LN_COBYLA;
opt.lower_bounds = 0*x0+xmin;
opt.upper_bounds = 0*x0+xmax;
opt.initial_step = 0*x0+xstep;
if ~isfield(opt,'xtol_rel'), opt.xtol_rel = 1e-8; end
opt.min_objective = @(x)feval(f,x(:));
nc = sum(~cellfun(@isnumeric,constr));
opt.fc = cell(1,nc);
rhs = zeros(1,nc);
k = 1;
for i=1:nc
  while isnumeric(constr{k}), k=k+1; end
  if k+1<=length(constr) && isnumeric(constr{k+1}), rhs(i)=constr{k+1}; end
  opt.fc{i} = eval(['@(x)feval(constr{' num2str(k) '},x(:))-rhs(' num2str(i) ')']);
  k = k+1;
end
[x, fopt, retcode] = nloptmin(opt, x0);
x = reshape(x,d,n);