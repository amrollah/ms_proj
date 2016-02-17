
s.loadframe(5918);
x = replace_isolated_nan(2*s.curseg.cc-1);
x(~s.mfi.sm) = nan;
x = cv.resize(x,.3);
x(~isnan(x)) = 2*(x(~isnan(x))>0)-1;

tic;
x1 = vmlMakeSurf(x);
toc;
figure(2);subplot(1,2,1);surf(x1);subplot(1,2,2);imagesc(sign(x1)+5*sign(x));
