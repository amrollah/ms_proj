function coeff1 = wlsq(X,Y,W,lambda,lb,ub,A,b,Aeq,beq)
if nargin<3, W = []; end
if nargin<4, lambda = []; end
if nargin<5, lb = []; end
if nargin<6, ub = []; end
if nargin<7, A = []; end
if nargin<8, b = []; end
if nargin<9, Aeq = []; end
if nargin<10, beq = []; end

[K,n] = size(X);
n0 = n;
onesK = ones(K,1);
assert(all(size(Y)==[K 1]));

if isempty(W), W = onesK; end
assert(all(size(W)==[K 1]));
W = W/mean(W);

jconst = [];%find(all(X==1),1);
if isempty(jconst)
  Gamma = eye(n);
  Delta = zeros(n,1);
  jkeep = true(1,n);
else
  xstd = std(X);
  xstd(jconst) = 1;
  jkeep = (xstd>=eps);
  X = X(:,jkeep);
  xstd = xstd(jkeep);
  xmean = mean(X);
  xmean(jconst) = 0;
  n = size(X,2);
  X = (X-xmean(onesK,:))./xstd(onesK,:);
  ymean = mean(Y);
  ystd = std(Y);
  Y = (Y-ymean)/ystd;
  jnotconst = setdiff(1:n,jconst);
  Gamma = diag(ystd./xstd);
  Gamma(jconst,jnotconst) = -ystd*xmean(jnotconst)./xstd(jnotconst);
  Delta = zeros(n,1);
  Delta(jconst) = ymean;
end

XW = X.*W(:,ones(1,n));
YW = Y.*W;
if isempty(lambda), lambda = zeros(n);
elseif min(size(lambda))==1, lambda = diag(lambda); end
if any(size(lambda)~=n), error('lambda has incorrect size'); end
if isempty(lb) && isempty(ub) && isempty(A) && isempty(b)
  %[Q,R] = qr(XW,0);
  %coeff = R\(YW'*Q)';
  coeff = (XW'*XW+lambda)\(XW'*YW);
else
  if ~isempty(A) || ~isempty(b)
    assert(size(b,2)==1);
    assert(all(size(A)==[size(b,1) n0]));
    A = A(:,jkeep);
  else
    A = zeros(0,n); b = zeros(0,1);
  end
  if ~isempty(Aeq) || ~isempty(beq)
    assert(size(beq,2)==1);
    assert(all(size(Aeq)==[size(beq,1) n0]));
    Aeq = Aeq(:,jkeep);
  else
    Aeq = zeros(0,n); beq = zeros(0,1);
  end
  if ~isempty(lb)
    assert(length(lb)==n0);
    lb = lb(jkeep);
  else
    lb = zeros(1,n)-inf; 
  end
  if ~isempty(ub)
    assert(length(ub)==n0);
    ub = ub(jkeep);
  else
    ub = zeros(1,n)+inf; 
  end
  Abounds = [-eye(n); eye(n)];
  bbounds = [-lb(:); ub(:)];
  Abounds(isinf(bbounds),:) = [];
  bbounds(isinf(bbounds)) = [];
  A = [Abounds; A];
  b = [bbounds; b];
  opt = [];
  opt.extractbounds = 1;
  opt.qfactor = 1;
  [status,coeff] = qlopt(XW'*XW+lambda,-YW'*XW,...
    A*Gamma,b-A*Delta,Aeq*Gamma,beq-Aeq*Delta,[],[],opt);
  assert(status==0);
end

coeff1 = zeros(n0,1);
coeff1(jkeep) = Gamma*coeff(:)+Delta;
