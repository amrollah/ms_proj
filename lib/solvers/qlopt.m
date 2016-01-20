function [status,x,fopt] = qlopt(Q,k,A,b,Aeq,beq,lb,ub,options)
% qlopt: mex interface to QP solver ql0001 (by Klaus Schittkowski)
%
% [status,x,fopt] = qlopt(Q,k,A,b,Aeq,beq,lb,ub,options)
%
% Options structure
%    options.maxabs           [double>0 (default Inf)]
%    options.boundstol        [double>0 (default 1e-8)]
%    options.extractbounds    [0 (off) or 1 (on) (default = 1, use to extract bounds from A, b, Aeq, Beq)]
%    options.qfactor          [double>=0 (default 2, use 2 for optimizing x'Qx instead of 0.5x'Qx)]
%
% please make sure to extract bounds if A or Aeq contain rows with a single
% nonzero element
%
%  status : 0 - optimal, 1 - infeasible, 2- max iter exceeded, 3- error

% written by Jan Poland, ABB.CHCRC.C1

if nargin<9
  options = [];
  if nargin < 8
    ub = [];
    if nargin < 7
      lb = [];
      if nargin < 6
        beq = [];
        if nargin < 5
          Aeq = [];
          if nargin < 4
            b = [];
            if nargin < 3
              A = [];
              if nargin < 2
                help qlopt;return
              end
            end
          end
        end
      end
    end
  end
end

A = [Aeq;-A];
b = [-beq;b];
neq = length(beq);
nx = length(k);
  
if isfield(options,'extractbounds'), extractbounds = options.extractbounds; else extractbounds = 1; end
if isfield(options,'qfactor'), qfactor = options.qfactor; else qfactor = 2; end
if isfield(options,'maxabs'), maxabs = options.maxabs; else maxabs=inf; end
  
if extractbounds
  if isempty(lb), lb=ones(1,nx)*(-maxabs); else lb=max(lb,-maxabs); end
  if isempty(ub), ub=ones(1,nx)*(maxabs); else ub=min(ub,maxabs); end
  lb = lb(:);
  ub = ub(:);
  
  nzA = sum(A~=0,2);
  ii = find(nzA==1);
  
  [ii1,jj,aa]=find(A(ii,:));
  ii2 = ii(ii1);
  vv = -b(ii2)./aa;
  jlb = find((ii2<=neq) | (aa>0));
  jub = find((ii2<=neq) | (aa<0));
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

  if isfield(options,'boundstol'), boundstol = options.boundstol; else boundstol=1e-8; end
  if any(lb>ub+boundstol)
    error('conflicting bounds!');
  end
  j = find(lb>ub);
  if ~isempty(j), lb(j)=ub(j); end
end

if isempty(lb), lb = zeros(nx,1)-maxabs; end
if isempty(ub), ub = zeros(nx,1)+maxabs; end

if qfactor~=1, Q=qfactor*Q; end
if issparse(Q), Q=full(Q); end
if issparse(k), k=full(k); end
if issparse(A), A=full(A); end
if issparse(b), b=full(b); end
if issparse(lb), lb=full(lb); end
if issparse(ub), ub=full(ub); end
if isempty(A), A=zeros(0,nx); b=zeros(0,1); end
[status,x,fopt] = ql0001mex(Q,k,A,b,neq,lb,ub);
