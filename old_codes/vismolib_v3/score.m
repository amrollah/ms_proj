function y = score(x,v,jj)
j = prod(sign(x(jj)-v),2)<0;
if sum(j)==0, y=0; 
else y = sum(abs(diff(x(jj(j,:)),[],2)))/sum(j); end
