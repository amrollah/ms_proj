function showsdf(sdf,x0)
if nargin<2, x0 = sign(sdf); end
clf;
p = 200+55*sign(sdf);
d1 = 48*(sign(sdf)>0&x0<0);
d2 = 80*(sign(sdf)<0&x0>0);
image(uint8(cat(3,p+d2,p-d1,p-d1)));
hold on; contour(sdf,20); hold off;
set(gca,'clim',[min(sdf(:)) max(sdf(:))]);
colorbar vert;
