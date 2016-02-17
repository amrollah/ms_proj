function x = vmlNormalizedRedChannel(curseg)
x = double(curseg.x(:,:,1));
x = x-mean(x(curseg.sm));
x = x/std(x(curseg.sm));
x(~curseg.sm) = NaN;
