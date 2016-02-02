function [seq1,m] = vmlMotSeg2(seq,tt,conf)

[Ny,Nx,T] = size(seq);
assert(length(tt)==T);
assert(tt(end)==0);
have_v = any(~isnan(seq),3);
%remove holes inside
have_v = ~(find_connected(~have_v,1,1) | find_connected(~have_v,Ny,1) | ...
  find_connected(~have_v,1,Nx) | find_connected(~have_v,Ny,Nx));

seq = reshape(seq,numel(seq)/T,T);
J = double(have_v); J(:) = cumsum(have_v(:)); J(~have_v) = NaN;

nx1 = nnz(have_v); nx = 2*nx1;

A = []; b = [];
Q = Q4sum(nx,nx1+J(:,1:end-1),nx1+J(:,2:end),sqrt(conf.wsmooth2_m),-sqrt(conf.wsmooth2_m)) ...
  + Q4sum(nx,nx1+J(1:end-1,:),nx1+J(2:end,:),sqrt(conf.wsmooth2_m),-sqrt(conf.wsmooth2_m));

k = zeros(nx,1);
% A = [A;A1];
% b = zeros(size(A,1),1) + conf.r2b_dmax;

% jj = J(have_v);
% A = [A;A4doublesum(nx,jj,nx1+jj,1,-tt(1))];
% b = [b;ones(2*nx1,1)];
  
%cost function for target
for t=1:T
  bb = ~isnan(seq(:,t)); jj = J(bb);
  w = 1./(1+abs(seq(bb,t)));
  k(jj) = k(jj)-conf.wtarget2*w.*seq(bb,t);
  k(nx1+jj) = k(nx1+jj)+conf.wtarget2*tt(t)*w.*seq(bb,t);
  Q = Q+Q4sum(nx,jj,nx1+jj,sqrt(w*conf.wtarget2),-tt(t)*sqrt(w*conf.wtarget2));
end

ub = [inf(1,nx1) inf(1,nx1)]; lb = -ub;

opts = [];
opts.extractbounds = 0;
opts.qfactor = 1;
[status,x]=ipoptqp(Q,k,A,b,[],[],lb,ub,opts);
assert(any(status==[0 1]),['Error: IpOpt returned status = ' num2str(status)]);

m = nan(size(have_v));
m(have_v) = -x(nx1+1:2*nx1);
seq1 = seq;
seq1(have_v,:) = bsxfun(@minus,x(1:nx1),bsxfun(@times,x(nx1+1:2*nx1),tt(:)'));
seq1 = reshape(seq1,[size(have_v) T]);

end

function Q = Q4sum(nx,j1,j2,c1,c2)
[jj,cc] = prepare_idz(j1,j2,c1,c2);
Q = sparse(nx,nx);
j = sub2ind(size(Q),jj(:,1),jj(:,1)); Q(j)=Q(j)+cc(:,1).^2;
j = sub2ind(size(Q),jj(:,2),jj(:,2)); Q(j)=Q(j)+cc(:,2).^2;
j = sub2ind(size(Q),jj(:,1),jj(:,2)); Q(j)=Q(j)+prod(cc,2);
j = sub2ind(size(Q),jj(:,2),jj(:,1)); Q(j)=Q(j)+prod(cc,2);
end

function A = A4sum(nx,j1,j2,c1,c2)
[jj,cc] = prepare_idz(j1,j2,c1,c2);
m = size(jj,1);
A = sparse(m,nx);
A(sub2ind(size(A),(1:m)',jj(:,1))) = cc(:,1);
A(sub2ind(size(A),(1:m)',jj(:,2))) = cc(:,2);
end

function A = A4doublesum(nx,j1,j2,c1,c2)
[jj,cc] = prepare_idz(j1,j2,c1,c2);
m = size(jj,1);
A = sparse(2*m,nx);
A(sub2ind(size(A),(1:m)',jj(:,1))) = cc(:,1);
A(sub2ind(size(A),(1:m)',jj(:,2))) = cc(:,2);
A(sub2ind(size(A),(m+1:2*m)',jj(:,1))) = -cc(:,1);
A(sub2ind(size(A),(m+1:2*m)',jj(:,2))) = -cc(:,2);
end

function [jj,cc] = prepare_idz(j1,j2,c1,c2)
if nargin<3, c1 = []; c2 = []; end
jj = [j1(:) j2(:)];
cc = [c1(:) c2(:)];
jdel = any(isnan(jj),2);
jj(jdel,:)=[];
if size(cc,1)>1, cc(jdel,:)=[]; end
end
