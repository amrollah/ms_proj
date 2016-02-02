function [xg,gd,xh,hd,hdd]=DerivativeCheck(f,params,xl,xu,bHess,n,d)

if nargin<4 | isempty(xu)
  if nargin<3, dim=1; else dim=xl; end;
  xl = zeros(dim,1);
  xu = ones(dim,1);
else
  dim = length(xl);
end;

if nargin<7, d=1e-8; end;
if nargin<6 | isempty(n), n=1; end;
if nargin<5, bHess=1; end;
xl=min(xu,xl+d);
xu=max(xl,xu-d);

ri = ones(1,dim);
gd = -1;
hd = -1;
hdd = -1;

for i=1:n
  x = xl+rand(dim,1).*max((xu-xl),0);
  if length(d)==1, D=eye(dim)*d; else D=diag(d); end;
  x = [x x(:,ri)+D];
  if bHess
    [y,g,H] = feval(f,x,params{:});
  else
    [y,g] = feval(f,x,params{:});
  end;
  
  g = reshape(g,[dim dim+1]);
  ng = reshape(y(2:end)-y(1),[dim 1])./d;
  graddiff = max(abs(g(:,1)-ng));
  if graddiff>gd, gd=graddiff; xg=x(:,1); end;
  
  if bHess
    H = reshape(H,[dim dim dim+1]);
    nH = g(:,2:end)-g(:,ri);
    if length(d)==1, nH = nH/d; else nH=nH./(d(:,ri)'); end;
    hessdiff = max(max(abs(H(:,:,1)-nH)));
    hessdiffdiag = max(abs(diag(H(:,:,1))-diag(nH)));
    if hessdiff>hd, hd=hessdiff; xh=x(:,1); end;
    if hessdiffdiag>hdd, hdd=hessdiffdiag; end;
  end;
    
  if rem(i,1000)==0 fprintf('.'); end;
end;

if (n==1)
  disp(['y = ' num2str(y(1))]);
  disp(['gradmax = ' num2str(max(abs(g(:,1))))]);
  if bHess
    disp(['hessmax = ' num2str(max(max(abs(H(:,:,1)))))]);  
  end;
else
  disp([num2str(n) ' points checked']);
end;

disp(' ');  
disp(['x = ' vec2str(xg)]);
disp(['graddiff = ' num2str(gd)]);
if bHess
  disp(' ');
  disp(['x = ' vec2str(xh)]);
  disp(['hessdiff = ' num2str(hd)]);
  disp(['hessdiff diag = ' num2str(hdd)]);
end;