function x = vmlColorify(x,mask,channels,offset)
mask(isnan(mask)) = false;
x = reshape(x,[numel(mask) 3]);
x(mask,channels) = max(0,min(255,x(mask,channels) + offset));
x = reshape(x,[size(mask) 3]);
