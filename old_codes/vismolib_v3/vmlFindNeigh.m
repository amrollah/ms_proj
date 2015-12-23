function neighInfo = vmlFindNeigh(sm,r)
if nargin<2, r=1; end
[dy,dx] = ndgrid(-floor(r):floor(r),0:floor(r));
d = [dy(:) dx(:)];
d(1:floor(r)+1,:) = [];
d(sum(d.^2,2)>r^2,:)=[];
[m,n] = size(sm);
neighInfo = [];
for i=1:size(d,1)
  ktop = max(0,-d(i,1));
  kbot = max(0,d(i,1));
  n1 = n-d(i,2);
  j = find([false(ktop,n1); sm(ktop+1:m-kbot,1:n1); false(kbot,n1)]);
  jneigh = logical(sm(j+d(i,1)+m*d(i,2)));
  neighInfo = [neighInfo;j(jneigh) j(jneigh)+d(i,1)+m*d(i,2)]; %#ok<AGROW>
end
