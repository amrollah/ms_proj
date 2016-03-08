function y = feval(m,x)
if size(x,1)~=size(m.x,1)
  error(['input must be in dimension ' num2str(size(m.x,1))])
end
[d,n] = size(m.x);
y = zeros(1,size(x,2));
for i=1:size(x,2)
  if isfield(m.kernel,'sig')
    K = exp(-0.5*(calcdist(m.x,x(:,i),1)/m.kernel.sig).^2);
  else
    K = (m.x'*x(:,i)/d+m.kernel.b).^m.kernel.n;
  end
  y(i) = m.a*K;
end
y=normalize(y,m.range(1),m.range(2));
