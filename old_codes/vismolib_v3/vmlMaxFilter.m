function y = vmlMaxFilter(x,ksize)
y = minmaxfilt(x,ksize+[0 0],'max','same');
