function [seq1,m] = vmlMotSeg(seq,tt,conf)

[Ny,Nx,T] = size(seq);
assert(length(tt)==T);
assert(tt(end)==0);
have_v = any(~isnan(seq),3);
r2 = bsxfun(@plus,((1:Ny)'-(1+Ny)/2).^2,((1:Nx)-(1+Nx)/2).^2);
have_v(r2<max(r2(have_v))) = true;

seq = reshape(seq,numel(seq)/T,T);
J = double(have_v); J(:) = cumsum(have_v(:)); J(~have_v) = NaN;
nslack1 = size(prepare_idz(J(:,1:end-1),J(:,2:end)),1);
nslack2 = size(prepare_idz(J(1:end-1,:),J(2:end,:)),1);

nx1 = nnz(have_v); nx = 2*nx1 + (nslack1+nslack2)*(conf.wsmooth1_m>0);

A = []; b = [];
Q = Q4sum(nx,J(:,1:end-1),J(:,2:end),sqrt(conf.wsmooth),-sqrt(conf.wsmooth));
Q = Q + Q4sum(nx,J(1:end-1,:),J(2:end,:),sqrt(conf.wsmooth),-sqrt(conf.wsmooth));
Q = Q + Q4sum(nx,nx1+J(:,1:end-1),nx1+J(:,2:end),sqrt(conf.wsmooth2_m),-sqrt(conf.wsmooth2_m));
Q = Q + Q4sum(nx,nx1+J(1:end-1,:),nx1+J(2:end,:),sqrt(conf.wsmooth2_m),-sqrt(conf.wsmooth2_m));

j = sub2ind(size(Q),nx1+(1:nx1),nx1+(1:nx1)); Q(j)=Q(j)+conf.wzero2_m;

k = zeros(nx,1);
% A = [A;A1];
% b = zeros(size(A,1),1) + conf.r2b_dmax;

jj = J(have_v);
A = [A;A4doublesum(nx,jj,nx1+jj,1,-tt(1))];
b = [b;ones(2*nx1,1)];

if conf.wsmooth1_m>0
  A1 = [A4doublesum(nx,nx1+J(:,1:end-1),nx1+J(:,2:end),1,-1);...
    A4doublesum(nx,nx1+J(1:end-1,:),nx1+J(2:end,:),1,-1)];
  A1(sub2ind(size(A1),1:2*(nslack1+nslack2),...
    [2*nx1+(1:nslack1) 2*nx1+(1:nslack1) 2*nx1+nslack1+(1:nslack2) 2*nx1+nslack1+(1:nslack2)])) = -1;
  A = [A;A1];
  b = [b;zeros(2*(nslack1+nslack2),1)];
  k(2*nx1+(1:(nslack1+nslack2))) = conf.wsmooth1_m;
end

if conf.band>0
  for t=1:T-1
    bb = seq(:,t)==1; jj = J(bb);
    A = [A;A4sum(nx,jj,nx1+jj,-1,tt(1))];
    b = [b;zeros(length(jj),1)+1-conf.band];
    bb = seq(:,t)==-1; jj = J(bb);
    A = [A;A4sum(nx,jj,nx1+jj,1,-tt(1))];
    b = [b;zeros(length(jj),1)-1+conf.band];
  end
end
  
%cost function for target
for t=1:T
  bb = ~isnan(seq(:,t)); jj = J(bb);
  k(jj) = k(jj)-(conf.wtarget1+conf.wtarget2)*seq(bb,t);
  k(nx1+jj) = k(nx1+jj)+(conf.wtarget1+conf.wtarget2)*tt(t)*seq(bb,t);
  Q = Q+Q4sum(nx,jj,nx1+jj,sqrt(conf.wtarget2),-tt(t)*sqrt(conf.wtarget2));
end

ub = [ones(1,nx1) ones(1,nx1)*2/abs(tt(1))]; lb = -ub;
if conf.wsmooth1_m>0
 lb = [lb zeros(1,nslack1+nslack2)];
 ub = [ub inf(1,nslack1+nslack2)];
end

if conf.band>0
  bb = seq(:,T)==1; lb(J(bb)) = 1-conf.band;
  bb = seq(:,T)==-1; ub(J(bb)) = -1+conf.band;
end

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
