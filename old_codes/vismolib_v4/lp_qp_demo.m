N = 50;
p = randn(1,N);
for i=2:N, p(i) = 0.8*p(i-1)+p(i); end
figure(1);clf;plot(p);

nx = N + N-1;
nc = 2*(N-1);

lb = zeros(1,nx);
ub = [ones(1,N) inf(1,N-1)];

w = 1;
c = [w*p ones(1,N-1)];

A = sparse(nc, nx);
j = sub2ind(size(A),1:N-1,1:N-1); A(j) = -1;
j = sub2ind(size(A),1:N-1,2:N); A(j) = 1;
j = sub2ind(size(A),1:N-1,N+(1:N-1)); A(j) = -1;
j = sub2ind(size(A),N-1+(1:N-1),1:N-1); A(j) = 1;
j = sub2ind(size(A),N-1+(1:N-1),2:N); A(j) = -1;
j = sub2ind(size(A),N-1+(1:N-1),N+(1:N-1)); A(j) = -1;

b = zeros(nc,1);

tic;
[x_lp,~,status] = clp([],c,A,b,[],[],lb,ub);
disp(['elapsed with LP : ' num2str(toc)])
disp(status);

figure(1);clf;plot(1:N,p,'b',1:N,1-x_lp(1:N),'-r.');

%now QP

nx = N;
nc = 0;

lb = zeros(1,nx);
ub = ones(1,nx);

c = w*p;
Q = sparse(nx,nx);
Q(1) = 1; Q(end) = 1;
j = sub2ind(size(Q),2:N-1,2:N-1); Q(j) = 2;
j = sub2ind(size(Q),1:N-1,2:N); Q(j) = -1;
j = sub2ind(size(Q),2:N,1:N-1); Q(j) = -1;

options = [];
options.qfactor = 2;
options.extractbounds = 0;

% tic;
% [x_qp_clp,~,status] = clp(Q,c',[1 zeros(1,nx-1)],1,[],[],lb,ub,options);
% disp(['elapsed with QP / CLP: ' num2str(toc)])
% disp(status);

tic;
[status,x_qp_ipopt] = ipoptqp(Q,c',[],[],[],[],lb,ub,options);
disp(['elapsed with QP / IpOpt: ' num2str(toc)])
disp(status);

% tic;
% [status,x_qp_ql] = qlopt(Q,c,[],[],[],[],lb,ub,options);
% disp(['elapsed with QP / QLOPT: ' num2str(toc)])
% disp(status);

figure(1);clf;plot(1:N,p,'b',1:N,1-x_lp(1:N),'-r.',1:N,1-x_qp_ipopt,'k',1:N,1-round(x_qp_ipopt),'ko');
