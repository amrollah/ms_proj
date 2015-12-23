function z = gausskern(sigma, x, y)

invsigma = 1./sigma(:);
[D, n] = size(x);
if nargin<3
  mu = mean(x,2);
  x = bsxfun(@times,bsxfun(@minus,x,mu),invsigma);
  y = x;
else
  [d, m] = size(y);
  if d~=D, error('dimensions do not match'); end
  mu = (m/(n+m))*mean(y,2) + (n/(n+m))*mean(x,2);
  x = bsxfun(@times,bsxfun(@minus,x,mu),invsigma); 
  y = bsxfun(@times,bsxfun(@minus,y,mu),invsigma);
end
z = exp(-.5*max(0,bsxfun(@plus,sum(x.*x,1)',bsxfun(@minus,sum(y.*y,1),2*x'*y))));
