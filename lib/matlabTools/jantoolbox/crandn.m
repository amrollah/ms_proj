function y = crandn(m,n,cf)
%y = crandn(m,n,cf)
%generates n series of correlated N(0,1) variables,
%where cov(y(t,i),y(t-1,i)) = cf

y = randn(m,n);
for i=2:m
  y(i,:)=cf*y(i-1,:)+sqrt(1-cf^2)*y(i,:);
end

