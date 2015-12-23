function v = SmoothThresMap(vtarget,w,conf)
[m,n] = size(vtarget);
nv = numel(vtarget);
nx = nv+1;
jj = reshape(1:nv,m,n);

Q = sparse(nx,nx);
Q(sub2ind(size(Q),1:nv,1:nv))= w(:);

j = sub2ind(size(Q),jj(:,1:end-1),jj(:,1:end-1)); Q(j)=Q(j)+conf.wsmooth;
j = sub2ind(size(Q),jj(:,2:end),jj(:,2:end)); Q(j)=Q(j)+conf.wsmooth;
j = sub2ind(size(Q),jj(:,1:end-1),jj(:,2:end)); Q(j)=Q(j)-conf.wsmooth;
j = sub2ind(size(Q),jj(:,2:end),jj(:,1:end-1)); Q(j)=Q(j)-conf.wsmooth;
j = sub2ind(size(Q),jj(1:end-1,:),jj(1:end-1,:)); Q(j)=Q(j)+conf.wsmooth;
j = sub2ind(size(Q),jj(2:end,:),jj(2:end,:)); Q(j)=Q(j)+conf.wsmooth;
j = sub2ind(size(Q),jj(1:end-1,:),jj(2:end,:)); Q(j)=Q(j)-conf.wsmooth;
j = sub2ind(size(Q),jj(2:end,:),jj(1:end-1,:)); Q(j)=Q(j)-conf.wsmooth;

k = [-vtarget(:).*w(:);0];
ncdiff = 2*((m-1)*n+m*(n-1));
A = sparse(ncdiff+2*nv,nx);
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

A(sub2ind(size(A),i+(1:nv),nx+zeros(1,nv))) = 1;
A(sub2ind(size(A),i+(1:nv),1:nv)) = -1;
i = i+nv;
A(sub2ind(size(A),i+(1:nv),nx+zeros(1,nv))) = -1;
A(sub2ind(size(A),i+(1:nv),1:nv)) = 1;

b = [zeros(ncdiff,1)+conf.max_local_variation;zeros(nv,1);zeros(nv,1)+conf.max_total_variation];

opts = [];
opts.extractbounds = 0;
opts.qfactor = 1;
[status,x]=ipoptqp(Q,k,A,b,[],[],zeros(1,nx)-8,zeros(1,nx)+8,opts);
assert(status==0,['Error: IpOpt returned status = ' num2str(status)]);
v = reshape(x(1:nv),m,n);
