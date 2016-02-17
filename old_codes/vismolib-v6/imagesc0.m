function imagesc0(x,z)
if nargin<2, z = zeros(size(x)); end
imagesc(x);colorbar vert;
hold on;
contour(x-z,[0 0],'k');
hold off;
set(gca,'clim',[min(x(:)) max(x(:))]);
