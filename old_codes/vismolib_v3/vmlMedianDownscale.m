function [X1,xx,yy] = vmlMedianDownscale(X,ksize)
X = vmlMedianBlur(X,ksize);
size1 = floor(size(X)/ksize);
yy = ceil(ksize/2)+floor((size(X,1)-size1(1)*ksize)/2)+(0:size1(1)-1)*ksize;
xx = ceil(ksize/2)+floor((size(X,2)-size1(2)*ksize)/2)+(0:size1(2)-1)*ksize;
X1 = X(yy,xx,:);
