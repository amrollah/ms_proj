function showgap(r2b,conf)
n = size(r2b,1);
mi = min(r2b(:));
ma = max(r2b(:));
xx = mi+(ma-mi)*(0:conf.findsep_N)/conf.findsep_N;
YY = zeros(n,conf.findsep_N+2);
YY(sub2ind(size(YY),(1:n)',1+ceil((r2b(:,1)-mi)/(ma-mi)*conf.findsep_N)))=1;
YY(sub2ind(size(YY),(1:n)',2+floor((r2b(:,2)-mi)/(ma-mi)*conf.findsep_N))) = -1;
yy = sum(max(0,cumsum(YY,2)),1);
yy = yy(1:end-1);
jrelax = yy>=max(yy)*conf.findsep_relax;
sep = find(jrelax,1);
plot(r2b(:,1),1:n,'b.',r2b(:,2),1:n,'r.',xx,yy,'k',xx(jrelax),yy(jrelax),'g.',xx(sep),yy(sep),'go');
