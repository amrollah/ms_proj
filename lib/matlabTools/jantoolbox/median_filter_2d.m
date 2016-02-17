function y = median_filter_2d(x,r,zero_is_nan,circ_mask,max_nan_frac)
if nargin<3, zero_is_nan = 0; end
if nargin<4, circ_mask = 0; end
if nargin<5, max_nan_frac = 0.5; end
if ~isa(x,'uint8') && ~isa(x,'uint16'), error('x must be uint8 or uint16'); end
y = x;
if size(x,1)>=size(x,2)
  for i=1:size(x,3)
    y(:,:,i) = median_filter_2d_core(y(:,:,i),double(r),...
      double(zero_is_nan),double(circ_mask),double(max_nan_frac));
  end
else
  for i=1:size(x,3)
    y(:,:,i) = median_filter_2d_core(y(:,:,i)',double(r),...
      double(zero_is_nan),double(circ_mask),double(max_nan_frac))';
  end
end
