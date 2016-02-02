function phat = gpwfit(x)
%phat = gpwfit(x) returns maximum likelihood estimates of the parameters
%of the GPW distribution given the data
%
%Jan Poland, ABB.CH-RD.C1, 2013-10-31

x = x(:);

if any(x<0), error('distribution only defined on non-negative data'); end

% initial guess is exponential
k0 = 0;
kw0 = 1;
sigma0 = max(eps,mean(x));

% initial guess is Weibull
% phat_wbl  = wblfit(x);
% k0 = 0;
% kw0 = phat_wbl(2);
% sigma0 = phat_wbl(1);

opt = optimset;
opt.MaxFunEvals = 400;
opt.MaxIter = 200;
opt.TolBnd = 1e-6;
opt.TolFun = 1e-6;
opt.TolX = 1e-6;

% Maximize the log-likelihood 
phat = fminsearch(@negloglik,[k0 log(kw0) log(sigma0)],opt,x);
phat(2:3) = exp(phat(2:3));

function y = negloglik(p, x)
k = p(1);
kw = exp(p(2));
sigma = exp(p(3));
n = numel(x);
z = x./sigma;
if abs(k) > eps
  if k > 0 || max(z.^kw) < -1/k
    y = n*log(sigma/kw) + (1-kw)*sum(log(z)) + (1+1/k) * sum(log1p(k*z.^kw));
  else
    y = inf;
  end
else % Weibull distribution for k close to 0
  y = n*log(sigma/kw) + (1-kw)*sum(log(z)) + sum(z.^kw);
end
