function y = branin(x);
if nargin<1 | isempty(x), x=2; end;
if length(x(:))==1
  if x>=0, y = [-pi+5 12.275; pi+5 2.275; 9.42478+5 2.475]'/15;
  else y = [0 0]'; end;
  return;
end;
a=1; b=5.1/(4*pi^2); c=5/pi; d=6; e=10; f=1/(8*pi);
x = x*15;
x(1,:)=x(1,:)-5;
y = a*(x(2,:)-b*x(1,:).^2+c*x(1,:)-d).^2+e*(1-f)*cos(x(1,:))+e; 
%xopt = [-pi+5 12.275; pi+5 2.275; 9.42478+5 2.475]'/15;