function y = vmlSobelX(x)
y = cv.Sobel(double(x),'XOrder',1,'YOrder',0);

