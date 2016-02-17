function mesh1(t,yy,xx)
seq = evalin('base','seq');
seq1 = evalin('base','seq1');
if nargin<3, yy=1:size(seq,1); xx=1:size(seq,2); end
seq = seq(yy,xx,t);
seq1 = seq1(yy,xx,t);
clf;
mesh(seq1);
hidden off;
z0 = get(gca,'zlim'); z0 = z0(1);
hold on;
contour(seq1,[0 0],'g','ContourZLevel',z0);
contour(seq,[0 0],'k','ContourZLevel',z0);
hold off;
set(gca,'clim',[min(seq1(:)) max(seq1(:))]);
