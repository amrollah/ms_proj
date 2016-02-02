function showgrad(x,m)
[gx,gy]=gradient(x);
clf;
imagesc(m);colorbar vert;
hold on;
contour(x,[0 0],'k');
quiver(gx.*m,gy.*m);
hold off;
set(gca,'clim',[min(m(:)) max(m(:))]);
