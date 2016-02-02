function y = fTanh2M(x);
y = fTanh2(x) + randn(1,size(x,2))*0.05; %*size(x,1)
