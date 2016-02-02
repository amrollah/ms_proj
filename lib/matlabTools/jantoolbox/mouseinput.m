function [xx,yy] = mouseinput(xx,yy)

if nargin<2
  xx = [];
  yy = [];
end;
  
figure('Name','Mouse input','DoubleBuffer','on');
axis([0 1 0 1]);
while 1
  [x,y]=ginput(1);
  if x<0 | x>1 | y<0 | y>1, break;end;
  xx = [xx x];yy = [yy y];
  plot(xx,yy,'.');
  axis([0 1 0 1]);
end;
close;
