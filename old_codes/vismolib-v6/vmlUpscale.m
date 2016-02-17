function y = vmlUpscale(x,s)
y = cv.resize(extrapolate2nan(x),s);