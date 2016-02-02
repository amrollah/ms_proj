function v = vmlSmoothThresMap(thres,have_v,conf)

% first find smooth map of starting values
N = conf.ngrid;
J = have_v; J(:) = cumsum(have_v(:)); J(~have_v) = NaN;
nx = nnz(have_v)+(conf.w_boundary>0);
[Q,A] = smooth_mats(nx,J(:,1:end-1),J(:,2:end),1);
[Q1,A1] = smooth_mats(nx,J(1:end-1,:),J(2:end,:),1);
Q = Q+Q1;
k = zeros(1,nx);
A = [A;A1];
b = zeros(size(A,1),1) + conf.r2b_dmax;

if conf.w_boundary>0
  jj = J(thres.isouter & ~isnan(J));
  n = length(jj);
  j = sub2ind(size(Q),jj,jj); Q(j)=Q(j)+conf.w_boundary/n;
  j = sub2ind(size(Q),nx,nx); Q(j)=Q(j)+conf.w_boundary;
  j = sub2ind(size(Q),jj,jj*0+nx); Q(j)=Q(j)-conf.w_boundary/n;
  j = sub2ind(size(Q),jj*0+nx,jj); Q(j)=Q(j)-conf.w_boundary/n;
end

%initial QP: try to hit local minima // NO LONGER USED
% Q1 = Q; k1 = k;
% w0 = 3;
jj = J(~isnan(thres.j0));
j = sub2ind(size(Q1),jj,jj); 
% Q1(j)=Q1(j)+w0; 
Q(j)=Q(j)+conf.w_local_min;
vv = thres.zzh(thres.j0(~isnan(thres.j0)));
% k1(jj) = -w0*vv; 
k(jj) = -conf.w_local_min*vv;

%defined prior on those without local minimum
jj = J(isnan(thres.j0)); jj(isnan(jj)) = [];
j = sub2ind(size(Q),jj,jj);
Q(j) = Q(j)+conf.w_quadprior;
k(jj) = k(jj)-conf.w_quadprior*conf.v_quadprior;

% minx = min(min(vv)); maxx = max(max(vv));
% lb = zeros(N)+minx; %lb(~isnan(thres.vlb)) = thres.vlb(~isnan(thres.vlb));
% ub = zeros(N)+maxx; %ub(~isnan(thres.vub)) = thres.vub(~isnan(thres.vub));
% lb = lb(logical(have_v)); ub = ub(logical(have_v));
% if conf.w_boundary>0, lb = [lb;minx]; ub = [ub;maxx]; end

opts = [];
opts.extractbounds = 0;
opts.qfactor = 1;
% [status,x0]=ipoptqp(Q1,k1,A,b,[],[],lb,ub,opts);
% assert(any(status==[0 1]),['Error: IpOpt returned status = ' num2str(status)]);

minx = thres.r2b_minmax(1); maxx = thres.r2b_minmax(2);
lb = zeros(N)+minx; %lb(~isnan(thres.vlb)) = thres.vlb(~isnan(thres.vlb));
ub = zeros(N)+maxx; %ub(~isnan(thres.vub)) = thres.vub(~isnan(thres.vub));
lb = lb(logical(have_v)); ub = ub(logical(have_v));
if conf.w_boundary>0, lb = [lb;minx]; ub = [ub;maxx]; end

has_tanh = squeeze(all(~isnan(thres.ww),1));
tanh_ix = repmat(J(has_tanh)',size(thres.ww,1),1);
tanh_w = thres.ww(:,has_tanh);
tanh_spec = [tanh_ix(:) repmat(1/conf.r2b_resolution,numel(tanh_ix(:)),1) ...
  repmat(-thres.zz(:)/conf.r2b_resolution,nnz(has_tanh),1) tanh_w(:)];

minxstart = min([conf.v_quadprior thres.zzh(thres.j0(~isnan(thres.j0)))]);
maxxstart = max([conf.v_quadprior thres.zzh(thres.j0(~isnan(thres.j0)))]);
n = ceil((maxxstart-minxstart)/conf.r2b_multistart_step)+1;
status = -1; fopt = inf;
for i=1:n
%   if i==0, x0_ = x0; else
  x0_ = zeros(nx,1)+minxstart+conf.r2b_multistart_step*(i-1)+(maxxstart-minxstart-(n-1)*conf.r2b_multistart_step)/2;
%   end
  [status1,x1,fopt1]=ipoptqp_tanh(x0_,tanh_spec,Q*conf.wsmooth,k*conf.wsmooth,A,b,[],[],lb,ub,opts);
  %disp([i x0_(1) status1 fopt1]);
  if any(status1==[0 1]) && fopt1<fopt
    status=status1; x=x1; fopt=fopt1;
  end
end

assert(any(status==[0 1]),['Error: IpOpt returned status = ' num2str(status)]);

v = have_v; v(logical(have_v)) = x(1:nnz(have_v)); v(~have_v) = NaN;

end

function [Q,A] = smooth_mats(nx,j1,j2,wsmooth)
jj = [j1(:) j2(:)];
jj(any(isnan(jj),2),:)=[];
Q = sparse(nx,nx);
j = sub2ind(size(Q),jj(:,1),jj(:,1)); Q(j)=Q(j)+wsmooth;
j = sub2ind(size(Q),jj(:,2),jj(:,2)); Q(j)=Q(j)+wsmooth;
j = sub2ind(size(Q),jj(:,1),jj(:,2)); Q(j)=Q(j)-wsmooth;
j = sub2ind(size(Q),jj(:,2),jj(:,1)); Q(j)=Q(j)-wsmooth;
m = size(jj,1);
A = sparse(2*m,nx);
A(sub2ind(size(A),(1:m)',jj(:,1))) = 1;
A(sub2ind(size(A),(1:m)',jj(:,2))) = -1;
A(sub2ind(size(A),(m+1:2*m)',jj(:,1))) = -1;
A(sub2ind(size(A),(m+1:2*m)',jj(:,2))) = 1;
end
