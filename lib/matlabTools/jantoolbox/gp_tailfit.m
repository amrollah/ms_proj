function [loc,p,xi,sigma] = gp_tailfit(x)
%[xi,sigma] = gp_fit(x,xi)
%fit Generalized Pareto distribution to data
%second argument xi is optional
%
%2014-01-08, Jan Poland, ABB.CH-RD.C1

N = 10;

x = sort(x(:),1,'descend');
j = round(length(x)*0.05);
% pp = ((j:-1:1)-.5)/j;
[xi,sigma] = gpfit1(x(1:j)-x(j+1));
% nll1 = [nll(x(1:N)-x(j+1),xi,sigma)/N nll(x(1:j)-x(j+1),xi,sigma)/j]-log(sigma);
% nlle = [nll_expected(pp(1:N),xi)/N nll_expected(pp,xi)/j];
% (nll1-nlle)./nlle
qe = q_expected(((j:-1:j-N+1)-.5)/j,xi)';
err = mean(abs((x(1:N)-x(j+1))/sigma-qe)./qe)
1;

function [xi,sigma] = gpfit1(x)
x = max(x,eps);
xis = [0 .1];
nlls = [nll(x,xis(1)) nll(x,xis(2))];
if nlls(1)<nlls(2), xi = fminbnd(@(xi) nll(x,xi),1/(1-max(x)/mean(x))+eps,xis(2));
else
  while 1
    xis = [xis xis(end)*2]; %#ok<*AGROW>
    nlls = [nlls nll(x,xis(end))];
    if nlls(end)>nlls(end-1), break; end
    xis(1)=[]; nlls(1)=[];
  end
  xi = fminbnd(@(xi) nll(x,xi),xis(1),xis(3));
end
sigma = get_sigma(x,xi);

function y=nll(x,xi,sigma)
if nargin<3, sigma = get_sigma(x,xi); end
if isnan(sigma), y = inf; return; end
z = x/sigma;
if abs(xi)>eps
  if xi > 0 || max(z)<-1/xi
    y = length(x)*log(sigma) + (1+1/xi) * sum(log1p(xi*z));
  else
    y = inf;
  end
else
  y = length(x)*log(sigma) + sum(z);
end

function sigma = get_sigma(x,xi)
j = max(1,floor(2*xi)+1);
B = fliplr(cumprod(1-2*xi./(length(x):-1:j)));
if isempty(B), 
  sigma=nan; else
  sigma = ((1+xi)*sum(B(2:end).*x(j+1:end)') + (j+1-xi)*B(1)*x(j))/sum(B);
end

function y = nll_expected(pp,xi)
if abs(xi) < eps
  y = sum(-log1p(-pp));
else
  y = (1+1/xi) * sum(log1p(expm1(-xi*log1p(-pp))));
end

function y = q_expected(pp,xi)
if abs(xi) < eps
  y = -log1p(-pp);
else
  y = expm1(-xi*log1p(-pp))/xi;
end


% function y=nll(x,xi)
% if abs(xi)>eps
%   sigma = mean(x)*(1-xi);
%   z = x/sigma;
%   if xi > 0 || max(z)<-1/xi
%     y = length(x)*log(sigma) + (1+1/xi) * sum(log1p(xi*z));
%   else
%     y = inf;
%   end
% else
%   sigma = mean(x);
%   y = length(x)*log(sigma) + sum(x/sigma);
% end
