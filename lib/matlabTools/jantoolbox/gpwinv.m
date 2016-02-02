function x = gpwinv(p,par)

k = par(1);
kw = par(2);
sigma = par(3);

pok = (0<p) & (p<1);
z = p+NaN;

if abs(k) < eps
  z(pok) = -log1p(-p(pok));
else
  z(pok) = expm1(-k*log1p(-p(pok)))/k;
end
if ~all(pok)
    % When k<0, the support is 0 <= (x-theta)/sigma <= -1/k
    z(p==0) = 0;
    if k<0, z(p==1)=-1/k; else z(p==1)=Inf; end
end

x = sigma*z.^(1/kw);
