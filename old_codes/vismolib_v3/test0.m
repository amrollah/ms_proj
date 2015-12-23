function test0(s,jj)
figure(111);clf;
for i=1:length(jj)
  j=jj(i);
  s.thresv(:,j)=nan;
  s.loadframe(0);
  subplot(3,length(jj),i);s.showCutoff(j);
  subplot(3,length(jj),i+length(jj));s.showR2BCutoff(j)
  subplot(3,length(jj),i+2*length(jj));s.showLumCutoff(j)
end