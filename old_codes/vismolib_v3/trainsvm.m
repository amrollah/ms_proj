function [w,b,btest]=trainsvm(s,tpred,gradb,train_frac,C)

if nargin<3, gradb = 0; end
if nargin<4, train_frac = 2/3; end
if nargin<5, C = 1000; end

ZZ = s.Z{s.getZidx(tpred)};
jj = find(all(~isnan(ZZ),1));
N = round(length(jj)*train_frac);

jtr = jj(1:N); 
Ztr = ZZ(:,jtr); 
Ptr = s.getP(s.ti(jtr) + tpred/86400); 
Ytr = s.cloudiness_label(s.ti(jtr) + tpred/86400);
jte = jj(N+1:end); 
Zte = ZZ(:,jte); 
Pte = s.getP(s.ti(jte) + tpred/86400);
Yte = s.cloudiness_label(s.ti(jte) + tpred/86400);
j = find(Ytr~=0);

d = size(Ztr,1);
N1 = sum(abs(Ytr));
if gradb(1)
  nx = d+2*N1;
  nc = N1+2*(N1-1);
else
  nx = d+1+N1;
  nc = N1;
end

Q = sparse(1:d,1:d,1,nx,nx);
k = [zeros(1,nx-N1) zeros(1,N1)+C];
if gradb(1)
  A1 = sparse([1:N1-1 1:N1-1],d+[1:N1-1 2:N1],[ones(1,N1-1) -ones(1,N1-1)],N1-1,nx);
  A = [-sparse(repmat(Ytr(j),d,1).*Ztr(:,j))' sparse(1:N1,1:N1,-Ytr(j))...
    -speye(N1);A1;-A1];
  b = [-ones(N1,1);repmat(diff(j(:))*gradb(1),2,1)];
else
  A = [-sparse([repmat(Ytr(j),d,1).*Ztr(:,j);Ytr(j)])' -speye(N1)];
  b = -ones(nc,1);
end
lb = [-inf(1,nx-N1) zeros(1,N1)];
ub = inf(1,nx);

opts = [];
opts.extractbounds = 0;
opts.qfactor = 1;

tic;
[status,x] = ipoptqp(Q,k,A,b,[],[],lb,ub,opts);
fprintf('SVM trained in %.1f s, status=%i\n',toc,status);

w = x(1:d)';
if gradb(1)
  b = interp1(j,x(d+1:d+N1),1:N);
  j1 = find(~isnan(b),1); b(1:j1-1) = b(j1);
  j1 = find(~isnan(b),1,'last'); b(j1+1:end) = b(j1);
else
  b = x(d+1);
end

Yte1 = Yte;
b1 = b(end);
aggressive = 1;
btest = Yte;
for i=1:length(Yte)
  btest(i) = b1;
  y = saturate(w*Zte(:,i)+b1);
  Yte1(i) = y;
  bgoal = b1;
  if aggressive || sign(Yte(i)*y)<0
    bgoal =  b1 + abs(Yte(i))*(Yte(i)-y); 
  end
  b1 = max(b1-gradb(end),min(b1+gradb(end),bgoal));
end
  
disp('training error');
disp_svm_err(Ytr,saturate(w*Ztr+b));
disp('test error');
disp_svm_err(Yte,Yte1);

