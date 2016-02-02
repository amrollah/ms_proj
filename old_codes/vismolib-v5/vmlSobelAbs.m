function y = vmlSobelAbs(x)
x = double(x);
y = max(abs(cv.Sobel(x,'XOrder',1,'YOrder',0)),...
        abs(cv.Sobel(x,'XOrder',0,'YOrder',1)));
