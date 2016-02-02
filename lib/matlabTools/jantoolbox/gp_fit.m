function [xi,sigma] = gp_fit(x,xi)
%[xi,sigma] = gp_fit(x,xi)
%fit Generalized Pareto distribution to data
%second argument xi is optional
%
%2014-01-08, Jan Poland, ABB.CH-RD.C1

x = sort(x(:),1,'descend');
if nargin<2 || isnan(xi),
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
end
sigma = get_sigma(x,xi);

function y=nll(x,xi)
sigma = get_sigma(x,xi);
if isnan(sigma), y = inf; else
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
  if isnan(y), y=inf; end
end

function sigma = get_sigma(x,xi)
j = max(1,floor(2*xi)+1);
B = fliplr(cumprod(1-2*xi./(length(x):-1:j)));
if isempty(B), sigma=NaN; else
  sigma = ((1+xi)*sum(B(2:end).*x(j+1:end)') + (j+1-xi)*B(1)*x(j))/sum(B);
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
