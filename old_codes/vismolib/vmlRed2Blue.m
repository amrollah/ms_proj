function r2b = vmlRed2Blue(x)
if ~isa(x,'double'), x = double(x)/255; end
if ndims(x)==2, x = permute(x,[1 3 2]); end
r2b = max(-8,min(8,log2(1e-8+x(:,:,1))-log2(1e-8+x(:,:,3))));
g2b = max(0,min(8,log2(1e-8+x(:,:,2))-log2(1e-8+x(:,:,3))));
j = r2b>0;
r2b(j) = min(r2b(j),3*g2b(j));

%y = (min(conf.r2b.max,max(conf.r2b.min,log2(red)-log2(1/255+blue)))-conf.r2b.min)/(conf.r2b.max-conf.r2b.min);
%y = 1+min(0,max(-1,log2(red+16/256)-log2(1/256+blue)));
