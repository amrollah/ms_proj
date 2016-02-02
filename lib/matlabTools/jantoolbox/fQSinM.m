function y = fQSinM(x);
y = fQSin(x) + randn(1,size(x,2))*0.05; %*size(x,1)
