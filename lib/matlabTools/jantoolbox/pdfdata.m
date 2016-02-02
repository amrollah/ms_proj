function p = pdfdata(xdata,xsupport,x)
xdata = sort(xdata(:));
nx = length(xdata);
xdata = [xsupport(1);xdata;xsupport(2)];
pxdata = [0 ((1:nx)-0.5)/nx 1]';
b = enforcemindistance(xdata,1e-7,1e-2);
xdata = xdata(b);
pxdata = pxdata(b);

%exponential fit to tails in case of infinite support
ptail = 0.001;
ntailmin = 10;
if ~isinf(xsupport(1)) 
  jleft=1; 
  muleft = 0; 
else
  [~,jleft]=min(abs(pxdata-ptail));
  jleft = max(jleft,1+ntailmin);
  muleft = mean(xdata(2:jleft-1)-xdata(jleft));
end
if ~isinf(xsupport(2)) 
  jright=length(xdata); 
  muright = 0; 
else
  [~,jright]=min(abs(pxdata-(1-ptail)));
  jright = min(jright,length(xdata)-ntailmin);
  muright = mean(xdata(jright+1:end-1)-xdata(jright));
end

p = zeros(size(x));
[x,jj]=sort(x(:));
j = findlastlower(x,xdata,1);
j(j==.5)=1; j=floor(j);
b = j<jleft; p(jj(b)) = pxdata(jleft)*exp(-(xdata(jleft)-x(b))/muleft)/max(muleft,1e-300); %left tail
b = j>jright; p(jj(b)) = (1-pxdata(jright))*exp(-(x(b)-xdata(jright))/muright)/max(muright,1e-300); %right tail
b = j>=jleft & j<=jright; j = j(b);
p(jj(b)) = (pxdata(j+1)-pxdata(j))./(xdata(j+1)-xdata(j));
