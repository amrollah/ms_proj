function showmot(m,dt,seq,seq_end)
if nargin<4, seq_end = seq(:,:,end); end
clf;
imagesc(m);colorbar vert;
hold on;
contour(seq(:,:,1),[0 0],'k');
contour(seq(:,:,end),[0 0],'r--');
contour(seq_end,[0 0],'r');
contour(seq_end+dt*m,[0 0],'g');
hold off;
set(gca,'clim',[min(m(:)) max(m(:))]);
