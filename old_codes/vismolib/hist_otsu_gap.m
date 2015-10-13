function hist_otsu_gap(res,hh)
if (~isempty(hh))
  ylim = max(hh(:,2));
  bar(hh(:,1),hh(:,2),1,'c');
  hold on;
  xx = [res([1 3],:);nan(1,2)];
  yy = [res(4,:);res(4,:)+res(2,:).*abs(diff(xx(1:2,:)));nan(1,2)];
  plot(xx(:),yy(:),'-r.');
  plot(res(end,1)+[0 0],[0 ylim],'-b.',res(5,1)+[0 0],[0 ylim],':b.',res(5,2)+[0 0],[0 ylim],':b.');
  hold off;
  %axis([get(gca,'xlim') 0 ylim]);
end
