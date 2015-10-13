function y = vmlRed2Blue(red,blue)
if ~isa(red,'double'), red = double(red)/255; end
if ~isa(blue,'double'), blue = double(blue)/255; end
y = max(-8,min(8,log2(1e-8+red)-log2(1e-8+blue)));
%y = (min(conf.r2b.max,max(conf.r2b.min,log2(red)-log2(1/255+blue)))-conf.r2b.min)/(conf.r2b.max-conf.r2b.min);
%y = 1+min(0,max(-1,log2(red+16/256)-log2(1/256+blue)));
