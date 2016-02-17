function y = vmlLocalVariation(x,ksize)
% y = vmlLocalVariation(x,ksize)
% Calculate local standard deviation of the difference of x and
% it's Gaussian blur with ksize.
% JP / ABB.CH-RD.C1 / 2016-01-21
x = x-vmlGaussianBlur(x,ksize);
y = std(x(~isnan(x)));
