function hist_with_gap(hh,cc)
if (~isempty(hh))
  bar(hh(1,:),hh(2,:),1,'c');
  hold on;
  xx1 = [cc(1,1) max(hh(1,:))];
  xx2 = [cc(1,2) min(hh(1,:))];
  plot(xx1,(xx1-xx1(1))*cc(2,1),'b',xx2,(xx2(1)-xx2)*cc(2,2),'b');
  if all(cc(2,3:4)>0)
    xx1 = [cc(1,4) max(hh(1,:))];
    xx2 = [cc(1,3) min(hh(1,:))];
    plot(xx1,(xx1-xx1(1))*cc(2,4),'r',xx2,(xx2(1)-xx2)*cc(2,3),'r');
  end
  axis([get(gca,'xlim') 0 max(hh(2,:))+1]);
  hold off;
end
