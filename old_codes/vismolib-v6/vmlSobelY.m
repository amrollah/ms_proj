function y = vmlSobelY(x)
y = cv.Sobel(double(x),'XOrder',0,'YOrder',1);
