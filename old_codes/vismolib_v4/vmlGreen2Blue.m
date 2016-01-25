function g2b = vmlGreen2Blue(x)
if ~isa(x,'double'), x = double(x)/255; end
if ndims(x)==2, x = permute(x,[1 3 2]); end
g2b = max(-8,min(8,log2(1e-8+x(:,:,2))-log2(1e-8+x(:,:,3))));
% g2b = log2(x(:,:,2))-log2(x(:,:,3));
