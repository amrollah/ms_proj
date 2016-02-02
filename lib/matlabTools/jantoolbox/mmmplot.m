function mmmplot(data,w,d)

if nargin<2, w=0.8; end;
if nargin<3, d=0.25; end;


m = mean(data);
s = std(data);
mi = min(data);
ma = max(data);
bar(m,w);hold on;
%set(gca,'XTickLabel','');
colormap([0.9 0.9 0.9]);
for i=1:length(m)
  plot(i(ones(size(data,1),1),1),data(:,i),'kx');
  gcline = line([i-0.1 i+0.1 i i i-0.1 i+0.1],[m(i)+s(i)...
      m(i)+s(i) m(i)+s(i) m(i)-s(i) m(i)-s(i) m(i)-s(i)]);
  set(gcline,'Color',[0.5 0.5 0.5],'LineWidth',1.5);
  gcline = line([i-d i+d i+d i-d i-d],...
    [mi(i) mi(i) ma(i) ma(i) mi(i)]);
  set(gcline,'Color',[0 0 0]);
end;
hold off;