function y = moving_avg_2d(x,r,circ_mask,max_nan_frac)
if nargin<3, circ_mask = 0; end
if nargin<4, max_nan_frac = 0.5; end
y = double(x);
if size(x,1)>=size(x,2)
  for i=1:size(x,3)
    y(:,:,i) = moving_avg_2d_core(y(:,:,i),...
      double(r),double(circ_mask),double(max_nan_frac));
  end
else
  for i=1:size(x,3)
    y(:,:,i) = moving_avg_2d_core(y(:,:,i)',...
      double(r),double(circ_mask),double(max_nan_frac))';
  end
end
