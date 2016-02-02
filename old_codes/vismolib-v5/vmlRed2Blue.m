function r2b = vmlRed2Blue(x)
if ~isa(x,'double'), x = double(x)/255; end
if ndims(x)==2, x = permute(x,[1 3 2]); end
r2b = max(-8,min(8,log2(1e-8+x(:,:,1))-log2(1e-8+x(:,:,3))));
