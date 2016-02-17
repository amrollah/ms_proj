function [seq1,m] = vmlMotSeg3(seq,tt,conf)

[Ny,Nx,T] = size(seq);
assert(length(tt)==T);
assert(tt(end)==0);
dt = diff(tt);
have_v = any(~isnan(seq),3);
r2 = bsxfun(@plus,((1:Ny)'-(1+Ny)/2).^2,((1:Nx)-(1+Nx)/2).^2);
have_v(r2<max(r2(have_v))) = true;

seq = reshape(seq,numel(seq)/T,T);
J = double(have_v); J(:) = cumsum(have_v(:)); J(~have_v) = NaN;
cumnslack = [0 cumsum(sum(~isnan(seq)))];
nx1 = nnz(have_v); nx = (T+1)*nx1+cumnslack(end);

A = []; b = [];

%smooth motion
Q = Q4sum(nx,J(:,1:end-1),J(:,2:end),sqrt(conf.wsmooth2_m),-sqrt(conf.wsmooth2_m));
Q = Q + Q4sum(nx,J(1:end-1,:),J(2:end,:),sqrt(conf.wsmooth2_m),-sqrt(conf.wsmooth2_m));

%smooth surfaces
for i=1:T
  Q = Q + Q4sum(nx,i*nx1+J(:,1:end-1),i*nx1+J(:,2:end),sqrt(conf.wsmooth2),-sqrt(conf.wsmooth2));
  Q = Q + Q4sum(nx,i*nx1+J(1:end-1,:),i*nx1+J(2:end,:),sqrt(conf.wsmooth2),-sqrt(conf.wsmooth2));
end

% fit of motion
jj = J(have_v);
j = sub2ind(size(Q),jj,jj); Q(j)=Q(j)+sum(dt.^2)*conf.wtarget2_m;
for i=1:T-1
  j = sub2ind(size(Q),i*nx1+jj,i*nx1+jj); Q(j)=Q(j)+conf.wtarget2_m;
  j = sub2ind(size(Q),(i+1)*nx1+jj,(i+1)*nx1+jj); Q(j)=Q(j)+conf.wtarget2_m;
  j = sub2ind(size(Q),(i+1)*nx1+jj,i*nx1+jj); Q(j)=Q(j)-conf.wtarget2_m;
  j = sub2ind(size(Q),i*nx1+jj,(i+1)*nx1+jj); Q(j)=Q(j)-conf.wtarget2_m;
  j = sub2ind(size(Q),(i+1)*nx1+jj,jj); Q(j)=Q(j)-dt(i)*conf.wtarget2_m;
  j = sub2ind(size(Q),jj,(i+1)*nx1+jj); Q(j)=Q(j)-dt(i)*conf.wtarget2_m;
  j = sub2ind(size(Q),i*nx1+jj,jj); Q(j)=Q(j)+dt(i)*conf.wtarget2_m;
  j = sub2ind(size(Q),jj,i*nx1+jj); Q(j)=Q(j)+dt(i)*conf.wtarget2_m;
end

%fit of target
k1 = seq(have_v,:); k1(isnan(k1)) = 0;
k = [zeros(1,nx1) -conf.wlo*k1(:)' conf.whi*ones(1,cumnslack(end))];
for i=1:T
  A1 = [A4doublesum(nx,i*nx1+J(:,1:end-1),i*nx1+J(:,2:end),1,-1);...
    A4doublesum(nx,i*nx1+J(1:end-1,:),i*nx1+J(2:end,:),1,-1)];
  A2 = [A4doublesum(nx,i*nx1+J(1:end-1,1:end-1),i*nx1+J(2:end,2:end),1,-1);...
    A4doublesum(nx,i*nx1+J(1:end-1,2:end),i*nx1+J(2:end,1:end-1),1,-1)];
  bb = ~isnan(seq(:,i)); jj = J(bb); nc1 = length(jj);
  A3 = sparse(nc1,nx);
  j = sub2ind(size(A3),1:nc1,i*nx1+jj'); A3(j) = -seq(bb,i);
  j = sub2ind(size(A3),1:nc1,(T+1)*nx1+cumnslack(i)+(1:nc1)); A3(j) = -1;
  A = [A;A1;A2;A3];
  b = [b;ones(size(A1,1),1);sqrt(2)*ones(size(A2,1),1);-conf.vsat*ones(nc1,1)];
end

lb = [-inf(1,nx1) -conf.ub*ones(1,T*nx1) zeros(1,cumnslack(end))];
ub = [inf(1,nx1) conf.ub*ones(1,T*nx1) inf(1,cumnslack(end))];

opts = [];
opts.extractbounds = 0;
opts.qfactor = 1;
[status,x]=ipoptqp(Q,k,A,b,[],[],lb,ub,opts);
assert(any(status==[0 1]),['Error: IpOpt returned status = ' num2str(status)]);

m = nan(size(have_v));
m(have_v) = x(1:nx1);
seq1 = seq;
seq1(have_v,:) = reshape(x(nx1+1:(T+1)*nx1),nx1,T);
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
