function x = dSimplex(d)
x = [-1 1];
D = 1;
for i=3:d  
  y = 1-1/(2*D^2);
  yy = sqrt(1-y^2);
  x = [yy*x zeros(i-2,1);-y+zeros(1,i-1) 1];
  D = yy*D;
end;