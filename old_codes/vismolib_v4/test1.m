function test1(s,jj)
figure(111);clf;
for i=1:length(jj)
  j=jj(i);
  s.thres.v(:,j)=nan;
  s.loadframe(0);
  subplot(3,length(jj),i);s.showR2B(j);
  subplot(3,length(jj),i+length(jj));s.showCutoff(j)
  subplot(3,length(jj),i+2*length(jj));s.showThresVmap(j);
end