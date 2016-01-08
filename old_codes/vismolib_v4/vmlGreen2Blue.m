function g2b = vmlGreen2Blue(x)
if ~isa(x,'double'), x = double(x)/255; end
if ndims(x)==2, x = permute(x,[1 3 2]); end
g2b = max(0,min(8,log2(1e-8+x(:,:,2))-log2(1e-8+x(:,:,3))));
