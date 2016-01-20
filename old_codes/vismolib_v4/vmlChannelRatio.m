function ch2ch = vmlChannelRatio(x,ch1,ch2)
if ~isa(x,'double'), x = double(x)/255; end
if nargin<3, ch2ch=x(:,:,ch1); return; end
if ndims(x)==2, x = permute(x,[1 3 2]); end
ch2ch = max(-8,min(8,log2(1e-8+x(:,:,ch1))-log2(1e-8+x(:,:,ch2))));
