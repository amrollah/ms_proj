function y = fGaussM(x);
y = fGauss(x) + randn(1,size(x,2))*0.03; %*size(x,1)
