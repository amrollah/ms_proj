function y = score1(x,xs,v,sm1)
x1 = x>v;
x1s = vmlSobelAbs(x1);x1s(~sm1)=0;x1s=vmlGaussianBlur(x1s,3);
y = xs(:)'*x1s(:)/sum(x1s(:));
