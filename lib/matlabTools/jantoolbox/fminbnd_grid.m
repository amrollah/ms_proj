function x = fminbnd_grid(f,xx)
yy = xx;
for i=1:length(xx), yy(i)=feval(f,xx(i)); end
[~,j]=min(yy);
if j==1, lb = xx(1); else lb=xx(j-1); end
if j==length(xx), ub = xx(end); else ub=xx(j+1); end
x = fminbnd(f,lb,ub);
