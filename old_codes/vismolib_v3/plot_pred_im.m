function plot_pred_im(s,j1,tpred)
n = ceil(sqrt(length(j1)));
clf;
for i = 1:length(j1)
  subplot(n,n,i);
  s.showPredIm(j1(i)+1,tpred);
  title(num2str(j1(i)+1));
end
