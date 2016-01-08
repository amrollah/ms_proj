function visualize_ngrids(x,n,xmin,xmax)
if nargin<3, xmin = nan; end
if nargin<4, xmax = nan; end
N = sqrt(length(x)/n);
if abs(N-round(N))>eps, error('incorrect data size'); end
N = round(N);
X1 = nan(N,n*N+(n-1)); X1(1:2,N+1) = [xmin;xmax];
for i=1:n, X1(1:N,(i-1)*(N+1)+(1:N)) = reshape(x((i-1)*N^2+(1:N^2)),N,N); end
imagesc(X1); axis off; colorbar vert;
mygray = [1 1 1]*.8;
for i=1:n-1,
  set(patch(i*(N+1)+[-.5 .5 .5 -.5 -.5],[0 0 N N 0]+.5,mygray),'EdgeColor',mygray);
end
