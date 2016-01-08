function y = vmlMinFilter(x,ksize)
y = minmaxfilt(x,ksize+[0 0],'min','same');
