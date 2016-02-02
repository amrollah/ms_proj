function p = fitskewnormal(x,prec)
if nargin<2, prec=0; end

n = length(x);
sskew = mean((x-mean(x)).^3)/var(x)^(3/2)/sqrt(6*(n-2)/(n+1)/(n+3));
%sample skewness is approximately normal distributed
if abs(sskew)<3
  p = [mean(x) std(x) 0];
  return;
end

y0 = nll(0,x); y = y0;
s00 = 0; s0 = 0; s = 1;
while 1
  y1 = nll(sign(sskew)*s,x);
  if y1>=y, break; end
  y = y1; s00 = s0; s0 = s; s = s*2;
end
s = fminbnd(@(a)nll(sign(sskew)*a,x),s00,s);
if prec>0, s = round(s/prec)*prec; end
[~,m,st]=nll(sign(sskew)*s,x);
p = [m st sign(sskew)*s];

% likelihood ratio test
% if gammainc(y0-y, 1/2)>0.99, p = [m st sign(sskew)*s];
% else p = [mean(x) std(x) 0]; end

function [y,m,s] = nll(a,x)
d=a/sqrt(1+a^2);
s = sqrt(var(x)/(1-2*d^2/pi));
m = mean(x)-s*d*sqrt(2/pi);
z = (x-m)/s;
y = length(x)*log(s)+sum(z.^2)/2-sum(log(erfc(-a*z/sqrt(2))));
