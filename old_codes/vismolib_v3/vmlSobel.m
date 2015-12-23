function y = vmlSobelY(x)
y = cv.Sobel(x,'XOrder',0,'YOrder',1);
