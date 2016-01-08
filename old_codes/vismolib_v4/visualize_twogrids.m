function visualize_twogrids(x,xmin,xmax)
if nargin<2, xmin = nan; end
if nargin<3, xmax = nan; end
N = sqrt(length(x)/2);
if abs(N-round(N))>eps, error('incorrect data size'); end
N = round(N);
X1 = [reshape(x(1:N^2),N,N) [xmin;xmax;nan(N-2,1)] reshape(x(N^2+1:end),N,N)];
imagesc(X1); axis off; colorbar vert;
mygray = [1 1 1]*.8;
set(patch(N+[.5 1.5 1.5 .5 .5],[0 0 N N 0]+.5,mygray),'EdgeColor',mygray);
