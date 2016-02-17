function [X1,xx,yy,ksize] = vmlMedianDownscale(X,s)
sz = size(X); sz = sz(1:2);
if length(s)==1 && s<=1, newsz = round(sz.*s); %scale given
else newsz = [s(1) s(end)];  %new size given
end
scale = sz./newsz;
ksize = round((scale-1)/2)*2+1;
if diff(ksize)~=0, warning('different ksize in x and y direction'); end
scale = min(scale);
ksize = max(ksize);
X = vmlMedianBlur(X,ksize);
pad = ((sz-1)-(newsz-1)*scale)/2;
yy = 1+pad(1)+(0:newsz(1)-1)*scale;
xx = 1+pad(2)+(0:newsz(2)-1)*scale;
X1 = X(round(yy),round(xx),:);
