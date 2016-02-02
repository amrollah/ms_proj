function qqdata(y,x,marker)
if nargin<3, marker='.'; end
y = sort(y);
x = sort(x);
py = ((1:length(y))-0.5)/length(y);
px = ((1:length(x))-0.5)/length(x);
xx = interp1(px,x,py,'linear','extrap');
mima = [min(xx(1),y(1)) max(xx(end),y(end))];
plot(xx,y,marker,mima,mima,'k--');
xlabel('quantiles of reference data');
ylabel('quantiles of examined data');
