function y = fWellenR(x);
y = fWellen(x) + randn(1,size(x,2))*0.15; %*size(x,1)
