function [y,f] = score2(r2b,lum,v,jj)
j = prod(sign(r2b(jj)-v),2)<0;
f = sum(j)/size(jj,1);
if sum(j)==0, y=0; 
else y = [sum(abs(diff(r2b(jj(j,:)),[],2))) sum(abs(diff(lum(jj(j,:)),[],2)))]/sum(j); end
