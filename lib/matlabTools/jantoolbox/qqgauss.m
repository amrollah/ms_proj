function qqgauss(x,mu,sigma,marker)
if nargin<4, marker='.'; end
x = sort(x(:));
if nargin<3
  mu = mean(x);
  sigma = std(x);
end
q1 = (-sqrt(2)*sigma).*erfcinv(2*((1:length(x))-0.5)/length(x)) + mu;
mima = [min(q1(1),x(1)) max(q1(end),x(end))];
plot(mima,mima,'r--',q1,x,marker);
xlabel('quantiles of normal distribution');
ylabel('quantiles of data');
