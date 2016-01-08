function [j1,j2] = vmlMakeIndexPair(n1,n2,nrep)
if nargin<3, nrep = 1; end
j1 = repmat(1:n1,n2,1);
j1 = j1(:);
j2 = repmat((1:n2)',1,n1);
j2 = j2(:);
if nrep>1
  j1 = bsxfun(@plus,j1,(0:nrep-1)*n1);
  j1 = j1(:);
  j2 = bsxfun(@plus,j2,(0:nrep-1)*n2);
  j2 = j2(:);
end
