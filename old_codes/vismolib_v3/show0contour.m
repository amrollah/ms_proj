function show0contour(s,i,thres)
if nargin<3, thres = 0; end
s.loadframe(i);
x = s.xcur.x; r2b = s.xcur.r2b;
x = vmlColorify(x,~s.mfi.sm,1:3,-255);
r2b(~s.mfi.sm) = nan;
image_level({x,r2b-thres},0,'g',1);
title([num2str(i) ' ' datestr(s.ti(i))]);
