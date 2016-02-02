function ax = phist(x,m)
if nargin<2, m=100; end
[nn,xx] = hist(x,m);
nn = nn/sum(nn);
[x1,n1] = stairs(xx,nn);
xx1 = (xx(1:end-1)+xx(2:end))/2;
dxx1 = xx1(2)-xx1(1);
ax = plotyy(x1,n1,[xx1(1)-dxx1 xx1 xx1(end)+dxx1],[0 cumsum(nn(1:end-1))/sum(nn) 1]);
set(ax,'ylim',[0 1]);
set(ax,'ytick',0:.1:1);
hold on; plot([xx(1) xx(end)],[0 0],'bd'); hold off;
set(ax,'fontsize',12);
grid on;
legend('histogram','min. / max.','cumulative distribution','location','best');
ylabel(ax(1),'relative #occurrences');
ylabel(ax(2),'cumulative distribution');

