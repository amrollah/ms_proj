function varargout = vmlSumsExceptDim(X,sz)
% varargout = vmlSumsExceptDim(X,sz)
% reshape X to sz and return all the sums of X along all dimensions except
% one
% JP / ABB.CH-RD.C1 / 2016-01-20

X = reshape(X,sz);
varargout = cell(1,length(sz));
for i = 1:length(sz)
  p = 1:length(sz); p(i) = [];
  varargout{i} = sum(reshape(permute(X,[i p]),sz(i),numel(X)/sz(i)),2);
end
