function res = jop(x,op,y)

if isempty(x) | isempty(y)
  res = [];
  return;
end;

sx = size(x);
sy = size(y);

if sx(1)>sy(1)
  if sy(1)~=1, error('Cannot adjust y!');end;
  ix1 = [1:sx(1)];
  iy1 = ones(1,sx(1));
elseif sx(1)<sy(1)
  if sx(1)~=1, error('Cannot adjust x!');end;
  ix1 = ones(1,sy(1));
  iy1 = [1:sy(1)];
else
  ix1 = [1:sx(1)];
  iy1 = ix1;
end;

if sx(2)>sy(2)
  if sy(2)~=1, error('Cannot adjust y!');end;
  ix2 = [1:sx(2)];
  iy2 = ones(1,sx(2));
elseif sx(2)<sy(2)
  if sx(2)~=1, error('Cannot adjust x!');end;
  ix2 = ones(1,sy(2));
  iy2 = [1:sy(2)];
else
  ix2 = [1:sx(2)];
  iy2 = ix2;
end;

res = feval(op,x(ix1,ix2),y(iy1,iy2));
