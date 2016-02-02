function X = squaredist(x,y);
if nargin<2, y=x; end;
[d,m] = size(x);
[d,n] = size(y);
Im = ones(1,m);
In = ones(1,n);
y = reshape(y,[d 1 n]);
X = reshape(sum((y(:,Im,:)-x(:,:,In)).^2,1),[m n]);
