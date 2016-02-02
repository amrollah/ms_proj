function y = ymess(x);
switch(size(x,1))
case 2, R=9.08;
case 3, R=14.89;
case 4, R=18.85;
case 5, R=21.91;
case 6, R=24.49;
case 7, R=26.74;
end;
y = fRadcosP(x);
y = y+randn(size(y))*0.01*R;
