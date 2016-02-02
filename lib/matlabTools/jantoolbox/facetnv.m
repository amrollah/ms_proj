function nvec = facetnv(x,k)
%calculate outer normal vectors to facets k of the convex
%hull of x
[nf, d] = size(k);
ri = ones(1,d-1);
nvec = zeros(d,nf);
for i=1:nf
  nv = null((x(:,k(i,2:end))-x(:,k(i,ri)))');
  sp = nv'*x-nv'*x(:,k(i,1));
  [v,j] = max(abs(sp));
  if sp(j)>0, nvec(:,i) = -nv; else nvec(:,i) = nv; end;  
end;