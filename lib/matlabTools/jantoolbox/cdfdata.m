function p = cdfdata(xdata,x)
xdata = sort(xdata(:));
pxdata = ((1:length(xdata))-0.5)/length(xdata);
p = min(1,max(0,interp1(xdata,pxdata,x,'linear','extrap')));
