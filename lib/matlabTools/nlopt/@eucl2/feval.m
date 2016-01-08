function y = feval(d,x)
if size(x,1)~=length(d.x)
  error(['input must be in dimension ' num2str(length(d.x))])
end
y = sqrt(sum((x-repmat(d.x,1,size(x,2))).^2,1));
