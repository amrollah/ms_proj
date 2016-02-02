function y = gp_quantile(xi,sigma,p)
if abs(xi)>eps, y = ((p^-xi)-1)*sigma/xi;
else y = -sigma*log(p); end
