function prel = vml_prel_clear_cavriglia(s,tt)
coeffs = [0.8535 0.2364 0.0227; 0.1342 0.1374 -0.2178];
c = coeffs(1,:)+coeffs(2,:)*cos(mod(s.data.day_of_year-1+11-183,365)/183*pi);
xx = s.getIrrClear(tt); xx = xx/max(xx);
prel = max(0.5*xx,1-c(1)*(1-xx)-c(2)*(1-xx).^2-c(3)*(1-xx).^3);
