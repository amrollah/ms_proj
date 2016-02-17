function x = vmlMakeSDF(cc,conf)

% wlo = .01;
% vsat = 0.5;
% wsmooth = 0.01; % not used (LP only)

x = 2*vmlDownscale(cc,conf.sz)-1;
have_v0 = ~isnan(x);

%remove holes inside
have_v = ~(find_connected(~have_v0,1,1) | find_connected(~have_v0,size(x,1),1) | ...
  find_connected(~have_v0,1,size(x,2)) | find_connected(~have_v0,size(x,1),size(x,2)));

x0 = x(have_v); x0 = x0(:);
J = double(have_v); J(:) = cumsum(have_v(:)); J(~have_v) = NaN;

nx1 = nnz(have_v); nxslack = nnz(have_v0); nx = nx1+nxslack;

% Q = Q4sum(nx,J(:,1:end-1),J(:,2:end),sqrt(wsmooth),-sqrt(wsmooth));
% Q = Q + Q4sum(nx,J(1:end-1,:),J(2:end,:),sqrt(wsmooth),-sqrt(wsmooth));

k = [zeros(1,nx1) ones(1,nxslack)];
k(J(have_v0)) = -conf.wlo*x0(J(have_v0));

A1 = [A4doublesum(nx,J(:,1:end-1),J(:,2:end),1,-1);...
    A4doublesum(nx,J(1:end-1,:),J(2:end,:),1,-1)];
A2 = [A4doublesum(nx,J(1:end-1,1:end-1),J(2:end,2:end),1,-1);...
    A4doublesum(nx,J(1:end-1,2:end),J(2:end,1:end-1),1,-1)];
A3 = [sparse(nxslack,nx1) -speye(nxslack)];
A3(sub2ind(size(A3),(1:nxslack)',J(have_v0))) = -x0(J(have_v0));
A = [A1;A2;A3];
b = [ones(size(A1,1),1);sqrt(2)*ones(size(A2,1),1);-conf.vsat*ones(nxslack,1)];

lb = [zeros(1,nx1)-conf.ub zeros(1,nxslack)];
ub = [zeros(1,nx1)+conf.ub inf(1,nxslack)];

[x1,~,status] = clp([],k,A,b,[],[],lb,ub,struct('solver',2));
%[status,x1] = ipoptqp(Q,k,A,b,[],[],lb,ub);

assert(status==0,'optimization failed');

x(have_v) = x1(1:nx1);

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

% function Q = Q4sum(nx,j1,j2,c1,c2)
% [jj,cc] = prepare_idz(j1,j2,c1,c2);
% Q = sparse(nx,nx);
% j = sub2ind(size(Q),jj(:,1),jj(:,1)); Q(j)=Q(j)+cc(:,1).^2;
% j = sub2ind(size(Q),jj(:,2),jj(:,2)); Q(j)=Q(j)+cc(:,2).^2;
% j = sub2ind(size(Q),jj(:,1),jj(:,2)); Q(j)=Q(j)+prod(cc,2);
% j = sub2ind(size(Q),jj(:,2),jj(:,1)); Q(j)=Q(j)+prod(cc,2);
% end
% 
% function A = A4sum(nx,j1,j2,c1,c2)
% [jj,cc] = prepare_idz(j1,j2,c1,c2);
% m = size(jj,1);
% A = sparse(m,nx);
% A(sub2ind(size(A),(1:m)',jj(:,1))) = cc(:,1);
% A(sub2ind(size(A),(1:m)',jj(:,2))) = cc(:,2);
% end
