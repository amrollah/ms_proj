function [status,x,fopt] = scip(c,A,b,Aeq,beq,lb,ub,vt,options,xinitial)
% SCIP Interface to MILP solver MEXSCIP
%
% [status,x,fopt] = scip(c,A,b,Aeq,beq,lb,ub,vt,options,xinitial)
%
% the first 8 arguments are mandatory
% c, A, b, Aeq, beq define the MILP in the standard form
%   min c'*x
%   s.t   A*x <  b, Aeq*x == beq, lb < x < ub
% lb and ub are lower and upper bounds on the x
% vt is a vector containing the variable types
%  0 - continuous (double)
%  1 - binary
%  2 - integer
%
% Options structure
%    options.maxnumseconds    [int>=0 (default 3600)]
%    options.verbose          [0|1|... (default 0)]
%    options.objgap           [double>0 (default 1e-6)]
%    options.maxabs           [double>0 (default Inf)]
%    options.boundstol        [double>0 (default 1e-8)]
%    options.extractbounds    [0 (off) or 1 (on) (default = 1, use to extract bounds from A, b, Aeq, Beq)]
%
% xinitial is presently not used
%
% please make sure to extract bounds if A or Aeq contain rows with a single
% nonzero element
%
% output
%  x      : solution
%  fopt   : optimal function value
%  status : 0 - optimal, 1 - infeasible, 2- unbounded

% Author of the original code: Johan Löfberg ETH Zürich.
% Extended by Jan Poland, ABB.CHCRC.C1

% **************************
% Check input
% **************************

if isstruct(c)
  if nargin<2, A = []; end
  if nargin<3, b = []; end  
  if isempty(A), options = []; else options = A; end
  xinitial = b;
  s = c;
  c=spmat(s.k); A=spmat(s.A); b=spmat(s.b); Aeq=spmat(s.Aeq); beq=spmat(s.beq); 
  vt=s.vt; lb=[]; ub=[];
else
  if nargin<10, xinitial = []; end
  if nargin<9, options = []; end
  if nargin<8, help scip; return; end
end

if isempty(c), c = ones(size(A,2),1); end

A = [Aeq;A];
b = [beq;b];
neq = length(beq);

nx = length(c);
  
if isfield(options,'maxabs'), maxabs = options.maxabs; else maxabs=inf; end
if isempty(lb), lb=ones(1,nx)*(-maxabs); else lb=max(lb,-maxabs); end
if isempty(ub), ub=ones(1,nx)*(maxabs); else ub=min(ub,maxabs); end
lb = lb(:);
ub = ub(:);

lb(vt==1) = max(lb(vt==1),0);
ub(vt==1) = min(ub(vt==1),1);

if ~isfield(options,'extractbounds') || options.extractbounds
  nzA = sum(A~=0,2);
  ii = find(nzA==1);
  
  [ii1,jj,aa]=find(A(ii,:));
  ii2 = ii(ii1);
  vv = b(ii2)./aa;
  jlb = find((ii2<=neq) | (aa<0));
  jub = find((ii2<=neq) | (aa>0));
  [vvlb,jjlb] = sort(vv(jlb));
  [vvub,jjub] = sort(vv(jub),1,'descend');
  j1lb = jj(jlb(jjlb));
  j1ub = jj(jub(jjub));
  lb(j1lb) = max(lb(j1lb), vvlb);
  ub(j1ub) = min(ub(j1ub), vvub);
  ii = find(nzA<2);
  A(ii,:)=[];
  b(ii)=[];
  neq = neq-sum(ii<=neq);
end

if isfield(options,'boundstol'), boundstol = options.boundstol; else boundstol=1e-8; end
if any(lb>ub+boundstol)
  error('conflicting bounds!');
end
j = find(lb>ub);
if ~isempty(j), lb(j)=ub(j); end

% **************************
% sparse format
% **************************
A = A';
cmatcnt = sum(A ~= 0,1);
cmatbeg = full(cumsum([0 cmatcnt]));
cmatbeg = cmatbeg(:)';
nzA = find(A);
cmatind = full(rem(nzA-1,size(A,1))');
cmatind = cmatind(:)';
cmatval = full(A(nzA));
cmatval = cmatval(:)';

c = full(c(:))';
b = full(b(:))';
lb = full(lb(:)');
ub = full(ub(:)');

[status,x,fopt] = mexscip(cmatbeg,cmatind,cmatval,c,b,neq,lb,ub,xinitial,vt,options);
