function f=rosen(x)
f = 2*sum((x(1:end-1,:).^2 - x(2:end,:)).^2,1) + ...
  sum((x(1:end-1,:)-1).^2,1);