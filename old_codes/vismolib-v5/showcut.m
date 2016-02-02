function showcut(seq,seq1,ii,jj)
figure(888);clf;hold on;
if length(ii)>length(jj), xx=ii; else xx=jj; end
ncc = size(seq,3);
cc = (0:ncc-1)'/(ncc-1)*pi;
set(gca,'ColorOrder',[max(0,cos(cc)) sin(cc) max(0,-cos(cc))]);
plot(xx,squeeze(seq(ii,jj,:)),':','linewidth',2);
% ax = gca;
% ax.ColorOrderIndex = 1;
plot(xx,squeeze(seq1(ii,jj,:)),'-');
hold off;