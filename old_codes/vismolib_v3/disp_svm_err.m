function disp_svm_err(Y,Y1)
j = isnan(Y) | isnan(Y1);
Y(j)=[]; Y1(j)=[];
N = sum(Y~=0);
Nclear = sum(Y>0);
Ncov = sum(Y<0);
Fcov = sum(Y>0 & Y1<0);
Fclear = sum(Y<0 & Y1>0);
F = sum(Y.*Y1<0);
fprintf('  class.: %i (%.1f%%) / f.clear: %i (%.1f%%) / f.cov: %i (%.1f%%) / N=%i\n',...
  F,100*F/N,Fclear,100*Fclear/Ncov,Fcov,100*Fcov/Nclear,N);
