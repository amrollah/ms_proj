function v = SmoothMapQP(vtarget,w,wsmooth,vlb,vub,dvub,plateaus)
if nargin<7, plateaus = zeros(0,4); end
[m,n] = size(vtarget);
nv = numel(vtarget);
jj = reshape(1:nv,m,n);

jslack = find(~isnan(plateaus(:,1)) & any(plateaus(:,3:4)>0,2));
jlo = find(plateaus(jslack,3)>0);
jhi = find(plateaus(jslack,4)>0);
nslack = length(jslack);
nx = nv+nslack;

Q = sparse(nx,nx);
Q(sub2ind(size(Q),1:nv,1:nv))= w(:);

j = sub2ind(size(Q),jj(:,1:end-1),jj(:,1:end-1)); Q(j)=Q(j)+wsmooth;
j = sub2ind(size(Q),jj(:,2:end),jj(:,2:end)); Q(j)=Q(j)+wsmooth;
j = sub2ind(size(Q),jj(:,1:end-1),jj(:,2:end)); Q(j)=Q(j)-wsmooth;
j = sub2ind(size(Q),jj(:,2:end),jj(:,1:end-1)); Q(j)=Q(j)-wsmooth;
j = sub2ind(size(Q),jj(1:end-1,:),jj(1:end-1,:)); Q(j)=Q(j)+wsmooth;
j = sub2ind(size(Q),jj(2:end,:),jj(2:end,:)); Q(j)=Q(j)+wsmooth;
j = sub2ind(size(Q),jj(1:end-1,:),jj(2:end,:)); Q(j)=Q(j)-wsmooth;
j = sub2ind(size(Q),jj(2:end,:),jj(1:end-1,:)); Q(j)=Q(j)-wsmooth;

k = [-vtarget(:).*w(:);ones(nslack,1)];
k(w==0)=0;
ncdiff = 2*((m-1)*n+m*(n-1));
A = sparse(ncdiff+length(jlo)+length(jhi),nx);
i = 0;
A(sub2ind(size(A),reshape(i+(1:m*(n-1)),m,n-1),jj(:,1:end-1))) = 1;
A(sub2ind(size(A),reshape(i+(1:m*(n-1)),m,n-1),jj(:,2:end))) = -1;
i = i+m*(n-1);
A(sub2ind(size(A),reshape(i+(1:m*(n-1)),m,n-1),jj(:,1:end-1))) = -1;
A(sub2ind(size(A),reshape(i+(1:m*(n-1)),m,n-1),jj(:,2:end))) = 1;
i = i+m*(n-1);
A(sub2ind(size(A),reshape(i+(1:(m-1)*n),m-1,n),jj(1:end-1,:))) = 1;
A(sub2ind(size(A),reshape(i+(1:(m-1)*n),m-1,n),jj(2:end,:))) = -1;
i = i+(m-1)*n;
A(sub2ind(size(A),reshape(i+(1:(m-1)*n),m-1,n),jj(1:end-1,:))) = -1;
A(sub2ind(size(A),reshape(i+(1:(m-1)*n),m-1,n),jj(2:end,:))) = 1;
i = i+(m-1)*n;

A(sub2ind(size(A),i+(1:length(jlo)),nv+jlo'))=-1;
A(sub2ind(size(A),i+(1:length(jlo)),jslack(jlo)'))=-plateaus(jslack(jlo),3);
i = i+length(jlo);
A(sub2ind(size(A),i+(1:length(jhi)),nv+jhi'))=-1;
A(sub2ind(size(A),i+(1:length(jhi)),jslack(jhi)'))=-plateaus(jslack(jhi),4);

b = [zeros(ncdiff,1)+dvub;...
  -plateaus(jslack(jlo),1).*plateaus(jslack(jlo),3);
  plateaus(jslack(jhi),2).*plateaus(jslack(jhi),4)];

opts = [];
opts.extractbounds = 0;
opts.qfactor = 1;
[status,x]=ipoptqp(Q,k,A,b,[],[],...
  [zeros(1,nv)+vlb zeros(1,nslack)],[zeros(1,nv)+vub zeros(1,nv)+inf],opts);
assert(status==0);
v = reshape(x(1:nv),m,n);
