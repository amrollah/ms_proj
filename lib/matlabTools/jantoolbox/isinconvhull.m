function [in,dist] = isinconvhull(xt,x,k,nv)
%#mex
%test if points xt are in the convex hull of x
%with facets k and outer normal vectors nv

d0 = sum(nv.*x(:,k(:,1)),1);
ri = zeros(1,size(k,1));
for i=1:size(xt,2)
  d = max(sum(nv.*xt(:,i+ri),1)-d0);
  if d<1e-7, in(i) = 1; else in(i) = 0; end;
  dist(i) = d;
end;
in = logical(in);
