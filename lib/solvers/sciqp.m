function [status,x,fopt] = sciqp(c,Q,A,QQ,blow,bupp,lb,ub,vt,options,xinitial)
% SCIP Interface to MIQCP solver MEXSCIP
%
% [status,x,fopt] = sciqp(c,Q,A,QQ,blow,bupp,lb,ub,vt,options,xinitial)
%
% the first 9 arguments are mandatory
%
% min c'*x + x'Qx
% s.t   blow <= x'*QQ*x + A*x <=  bupp, lb < x < ub
%
% QQ is represented as a cell array of matrices, such that each
% constraint reads blow(i) <= x'*QQ{i}*x + A(i,:)*x <= bupp(i)
% For linear constraints, the corresponding QQ{i} can be empty
%
% Options structure
%    options.maxnumseconds    [int>=0 (default 3600)]
%    options.verbose          [0|1|... (default 0)]
%    options.objgap           [double>0 (default 1e-6)]
%    options.maxabs           [double>0 (default Inf)]
%    options.boundstol        [double>0 (default 1e-8)]
%    options.extractbounds    [0 (off) or 1 (on) (default = 1, use to extract bounds from A, b, Aeq, Beq)]
%
% please make sure to extract bounds if A or Aeq contain rows with a single
% nonzero element
%
% output
%  x      : solution
%  fopt   : optimal function value
%  status : scip return code


if isstruct(c)
  if nargin<2, Q = []; end
  if nargin<3, A = []; end  
  if isempty(Q), options = []; else options = Q; end
  if isfield(options,'maxabs'), maxabs = options.maxabs; else maxabs=inf; end
  xinitial = A;
  s = c;
  c=spmat(s.k); Q=spmat(s.Q); 
  A=[spmat(s.Aeq); spmat(s.A)]; 
  blow=[spmat(s.beq);zeros(s.b(1),1)-maxabs];
  bupp=[spmat(s.beq);spmat(s.b)];
  vt=s.vt; lb=[]; ub=[];
  QQ = {};
else
  if nargin<11, xinitial = []; end
  if nargin<10, options = []; end
  if nargin<9, help sciqp; return; end
end

nx = length(c);
  
if isfield(options,'maxabs'), maxabs = options.maxabs; else maxabs=inf; end
if isempty(lb), lb=ones(1,nx)*(-maxabs); else lb=max(lb,-maxabs); end
if isempty(ub), ub=ones(1,nx)*(maxabs); else ub=min(ub,maxabs); end
lb = lb(:);
ub = ub(:);

lb(vt==1) = max(lb(vt==1),0);
ub(vt==1) = min(ub(vt==1),1);

if ~iscell(QQ), QQ={QQ}; end
QQ = QQ(:)';
if length(QQ)<size(A,1), QQ = [QQ cell(1,size(A,1)-length(QQ))]; end
if length(QQ)>size(A,1), QQ = QQ(1:size(A,1)); end

if ~isempty(A) && (~isfield(options,'extractbounds') || options.extractbounds)
  nzA = sum(A~=0,2);
  ii = find((nzA==1)' & cellfun(@isempty,QQ));
  [ii1,jj,aa]=find(A(ii,:));
  ii2 = ii(ii1);
  vvl = blow(ii2)./aa;
  vvu = bupp(ii2)./aa;
  vvl1 = vvl; 
  vvl(aa<0) = vvu(aa<0); vvu(aa<0) =vvl1(aa<0);
  [vvl,jjlb] = sort(vvl);
  [vvu,jjub] = sort(vvu,1,'descend');
  j1lb = jj(jjlb);
  j1ub = jj(jjub);
  lb(j1lb) = max(lb(j1lb), vvl);
  ub(j1ub) = min(ub(j1ub), vvu);
  ii = find(nzA<2);
  A(ii,:)=[];
  blow(ii)=[];
  bupp(ii)=[];
  QQ(ii)=[];
end

if isfield(options,'boundstol'), boundstol = options.boundstol; else boundstol=1e-8; end
if any(lb>ub+boundstol), error('conflicting bounds!'); end
j = find(lb>ub);
if ~isempty(j), lb(j)=ub(j); end
if any(blow>bupp+boundstol), error('conflicting bounds (blow>bupp)!'); end
j = find(blow>bupp);
if ~isempty(j), blow(j)=bupp(j); end

hasQ = 0;
if ~isempty(Q) && any(Q(:)~=0)
  hasQ = 1;
  A = [[A sparse(size(A,1),1)];c -1];
  blow = [blow; 0];
  bupp = [bupp; 0];
  QQ = [QQ {Q}];
  c = [sparse(1,nx) 1];
  lb = [lb;-maxabs];
  ub = [ub;maxabs];
end

A = A';
cmatcnt = sum(A ~= 0,1);
cmatbeg = full(cumsum([0 cmatcnt]));
cmatbeg = cmatbeg(:)';
nzA = find(A);
cmatind = full(rem(nzA-1,size(A,1))');
cmatind = cmatind(:)';
cmatval = full(A(nzA));
cmatval = cmatval(:)';

nnzQ = 0;
for i=1:length(QQ)
  QQ{i} = tril(QQ{i})+triu(QQ{i},1)';
  nnzQ = nnzQ+nnz(QQ{i});
end
cmatbegq = zeros(1,length(QQ)+1);
cmatindq1 = zeros(1,nnzQ);
cmatindq2 = zeros(1,nnzQ);
cmatvalq = zeros(1,nnzQ);
nnzQ = 0;
for i=1:length(QQ)
  [ii,jj,vv]=find(QQ{i});
  nnz1 = length(ii);
  cmatbegq(i+1) = cmatbegq(i)+nnz1;
  cmatindq1(nnzQ+(1:nnz1)) = ii-1;
  cmatindq2(nnzQ+(1:nnz1)) = jj-1;
  cmatvalq(nnzQ+(1:nnz1)) = vv;
  nnzQ = nnzQ+nnz1;
end

c = full(c(:))';
blow = full(blow(:))';
bupp = full(bupp(:))';
lb = full(lb(:)');
ub = full(ub(:)');

[status,x,fopt] = mexsciqp(cmatbeg,cmatind,cmatval,cmatbegq,cmatindq1,cmatindq2,cmatvalq,...
  c,blow,bupp,lb,ub,xinitial,vt,options);

if hasQ, x=x(1:end-1); end
