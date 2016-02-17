function y = f_pm(ccc,show)
if nargin<2, show = 0; end
ccc = reshape(ccc,2,3);
doy = evalin('base','doy');
xx1 = evalin('base','xx1');
yy1 = evalin('base','yy1');
y = 0;
ii = true(1,length(doy));
%ii([5 6 8]) = false;
for i=1:length(doy)
  c2 = ccc(1,:)+ccc(2,:)*cos(mod(doy(i)-1+11-183,365)/183*pi);
  xx = xx1{i}; yy = yy1{i};
  yp = max(0.5*xx,1-c2(1)*(1-xx)-c2(2)*(1-xx).^2-c2(3)*(1-xx).^3);
  y1 = mean((yp-yy).^2);
  if show 
    disp([num2str(i) ' ' num2str(y1)]); 
    if show>=2, figure(44+i);clf;plot(xx,yy,'b.',xx,yp,'r');end
  end
  if ii(i), y = y+y1; end
end
disp(y);
