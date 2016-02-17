function v = vmlSmoothThresMap(thres,have_v,conf)

% first find smooth map of starting values
N = conf.ngrid; N2=N^2;
J = have_v; J(:) = cumsum(have_v(:)); J(~have_v) = NaN;
nx1 = nnz(have_v); nx = nx1+(conf.w_boundary>0);
sz = size(J); J1 = reshape(1:prod(sz),sz);
[Q,A,tanh_pairs_spec] = neigh_info(1:sz(1),1:sz(2)-1,1:sz(1),2:sz(2),1);
[Q1,A1,tanh_pairs_spec1] = neigh_info(1:sz(1)-1,1:sz(2),2:sz(1),1:sz(2),1);
Q = Q + Q1; A = [A; A1]; 
tanh_pairs_spec = [tanh_pairs_spec; tanh_pairs_spec1];
b = zeros(size(A,1),1) + conf.r2b_dmax;
k = zeros(1,nx);
[Q1,A1] = neigh_info(1:sz(1)-1,1:sz(2)-1,2:sz(1),2:sz(2),.5);
Q = Q + Q1; A = [A; A1]; 
[Q1,A1] = neigh_info(1:sz(1)-1,2:sz(2),2:sz(1),1:sz(2)-1,.5);
Q = Q + Q1; A = [A; A1]; 
b = [b; zeros(size(A,1)-length(b),1) + conf.r2b_dmax*sqrt(2)];

if conf.w_boundary>0
  jj = J(thres.isouter & ~isnan(J));
  n = length(jj);
  j = sub2ind(size(Q),jj,jj); Q(j)=Q(j)+conf.w_boundary/n;
  j = sub2ind(size(Q),nx,nx); Q(j)=Q(j)+conf.w_boundary;
  j = sub2ind(size(Q),jj,jj*0+nx); Q(j)=Q(j)-conf.w_boundary/n;
  j = sub2ind(size(Q),jj*0+nx,jj); Q(j)=Q(j)-conf.w_boundary/n;
end

% j_j0 = J(~isnan(thres.j0));
v_j0 = thres.zzh(thres.j0(~isnan(thres.j0)));

% %target on identified local minima
% j = sub2ind(size(Q1),j_j0,j_j0);
% wlm = vmlGaussianBlur(double(~isnan(thres.j0)),3);
% wlm = (1-1/9*max(0,min(1,(wlm(~isnan(thres.j0))-0.25)/0.75)))*conf.w_local_min;
% Q(j)=Q(j)+wlm;
% k(j_j0) = -wlm.*v_j0(:);
% 
% %initial QP: try to hit local minima // NO LONGER USED
% % Q1 = Q; k1 = k;
% % w0 = 3;
% % Q1(j)=Q1(j)+w0; 
% % k1(j_j0) = -w0*v_j0; 
% 
% %prior on those without local minimum
% jj = J(isnan(thres.j0)); jj(isnan(jj)) = [];
% j = sub2ind(size(Q),jj,jj);
% Q(j) = Q(j)+conf.w_quadprior;
% k(jj) = k(jj)-conf.w_quadprior*conf.v_quadprior;

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

status = -1; fopt = inf;

minxstart = min([conf.prior_thres(2:3) v_j0]);
maxxstart = max([conf.prior_thres(2:3) v_j0]);
n = ceil((maxxstart-minxstart)/conf.r2b_multistart_step)+1;
x0 = thres.j0;
x0(~isnan(x0)) = v_j0;
while any(isnan(x0(:))), x0=extrapolate2nan(x0); end
x0 = x0(logical(have_v));
if conf.w_boundary>0, x0 = [x0;mean(x0)]; end
% sigma = (maxxstart-minxstart)/6;
for i=0:n
  if i==0, x0_ = x0; else
    x0_ = zeros(nx,1)+minxstart+conf.r2b_multistart_step*(i-1)+(maxxstart-minxstart-(n-1)*conf.r2b_multistart_step)/2;
  end
  [status1,x1,fopt1]=ipoptqp_tanh(x0_,tanh_spec,tanh_pairs_spec,Q*conf.wsmooth,k*conf.wsmooth,A,b,[],[],lb,ub,opts);
  if any(status1==[0 1]) && fopt1<fopt
    status=status1; x=x1; fopt=fopt1; x0=x;
  end
end
% for i=1:0
%   x0_ = max(minxstart,min(maxxstart,x0+randn(nx,1)*sigma*5/(i+5)));
%   [status1,x1,fopt1]=ipoptqp_tanh(x0_,tanh_spec,tanh_pairs_spec,Q*conf.wsmooth,k*conf.wsmooth,A,b,[],[],lb,ub,opts);
%   if any(status1==[0 1]) && fopt1<fopt
%     status=status1; x=x1; fopt=fopt1; x0=x;
%   end
% end

v = have_v; v(logical(have_v)) = x(1:nx1); v(~have_v) = NaN;
assert(any(status==[0 1]),['Error: IpOpt returned status = ' num2str(status1)]);


  function [Q,A,tanh_pairs_spec] = neigh_info(ry1,rx1,ry2,rx2,wsmooth)
    j1 = J(ry1,rx1); j2 = J(ry2,rx2);
    j11 = J1(ry1,rx1); j12 = J1(ry2,rx2);
    jj_ = [j1(:) j2(:) j11(:) j12(:)];
    jj_(any(isnan(jj_),2),:) = [];
    Q = sparse(nx,nx);
    j_ = sub2ind(size(Q),jj_(:,1),jj_(:,1)); Q(j_)=Q(j_)+wsmooth;
    j_ = sub2ind(size(Q),jj_(:,2),jj_(:,2)); Q(j_)=Q(j_)+wsmooth;
    j_ = sub2ind(size(Q),jj_(:,1),jj_(:,2)); Q(j_)=Q(j_)-wsmooth;
    j_ = sub2ind(size(Q),jj_(:,2),jj_(:,1)); Q(j_)=Q(j_)-wsmooth;
    m_ = size(jj_,1);
    A = sparse(2*m_,nx);
    A(sub2ind(size(A),(1:m_)',jj_(:,1))) = 1;
    A(sub2ind(size(A),(1:m_)',jj_(:,2))) = -1;
    A(sub2ind(size(A),(m_+1:2*m_)',jj_(:,1))) = -1;
    A(sub2ind(size(A),(m_+1:2*m_)',jj_(:,2))) = 1;
    ix_ = jj_(:,1:2)';
    s_ = [thres.p_tanh_cumsum(jj_(:,3)) thres.p_tanh_cumsum(jj_(:,4))]';
    b_ = [thres.p_tanh_cumsum(N2+jj_(:,3)) thres.p_tanh_cumsum(N2+jj_(:,4))]';
    w_ = [thres.p_tanh_cumsum(2*N2+jj_(:,3)) thres.p_tanh_cumsum(2*N2+jj_(:,4))]';
    wsp_ = min(thres.single_peak_crit(jj_(:,3)),thres.single_peak_crit(jj_(:,4))).*...
      max(thres.hist_left_open_crit(jj_(:,3)).*thres.hist_left_open_crit(jj_(:,4)),...
      thres.hist_right_open_crit(jj_(:,3)).*thres.hist_right_open_crit(jj_(:,4)));
    wsp_ = conf.scale_single_peak_crit*[wsp_ -wsp_]';
    tanh_pairs_spec = [ix_(:) s_(:) b_(:) w_(:).*wsp_(:)];
    tanh_pairs_spec(tanh_pairs_spec(:,4)==0,:) = [];
  end
end
