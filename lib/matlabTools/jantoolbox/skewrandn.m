function y = skewrandn(skew,m,n)
%y = skewrandn(skew,m,n)
%generates skew normal random variables
%default value for m and n is 1

if nargin<2, m=1; end;
if nargin<3, n=m; end;
y = randn(m,n);
j = (randn(m,n)>skew.*y);
y(j) = -y(j);
