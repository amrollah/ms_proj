function x = vmlColorify(x,mask,channels,v,overwrite)
if nargin<5, overwrite = 0; end
mask(isnan(mask)) = false;
x = reshape(x,[numel(mask) 3]);
if overwrite
  x(mask,channels) = v;
else
  x(mask,channels) = max(0,min(255,x(mask,channels) + v));
end
x = reshape(x,[size(mask) 3]);
