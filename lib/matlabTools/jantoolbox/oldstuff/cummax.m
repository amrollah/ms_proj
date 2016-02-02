function A = cummax(A,dim)
if nargin<2, dim=find(size(A)>1,1); end
if isempty(dim), return; end
if dim~=1
  perm=1:length(size(A)); 
  perm(dim)=1; perm(1)=dim;
  A = permute(A,perm);
end
for i=2:size(A,1), A(i,:)=max(A(i,:),A(i-1,:)); end
if dim~=1, A = permute(A,perm); end
