function [r,varargout]=decodeMessParams(nr,varargin)
nr = nr-1;
for i=1:min(nargin-1,nargout-1)
  a = varargin{i};
  k = length(a);
  varargout(i) = {a(rem(nr,k)+1)};
  nr = floor(nr/k);
end;
r = (nr==0);
