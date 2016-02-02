function drawGitter(x,y,z)

if nargin<3
  z = x;
  x = (1:(size(z,1)))';
  y = (1:(size(z,2)))';
else
  if min(size(x))==1
    x = x(:);
    y = y(:);
  elseif abs(x(end,1)-x(1,1))>abs(x(1,end)-x(1,1))
    x = x(:,1);
    y = y(1,:)';
  else
    x = x(1,:)';
    y = y(:,1);
  end;
end;

xx = [];
yy = [];
zz = [];

rix = ones(1,length(x));
riy = ones(1,length(y));

for i=1:length(x)
  xx = [xx x(i,riy) NaN];
  yy = [yy y' NaN];
  zz = [zz z(i,:) NaN];
end;
for i=1:length(y)
  xx = [xx x' NaN];
  yy = [yy y(i,rix) NaN];
  zz = [zz z(:,i)' NaN];
end;

plot3(xx,yy,zz,'k-');
