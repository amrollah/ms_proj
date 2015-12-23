function y = saturate(x)
y = x;
y(x<-1) = -1;
y(x>1) = 1;
