function [status,x,fopt] = ipoptqp_tanh(x0,tanh_spec,tanh_pairs_spec,Q,k,A,b,Aeq,beq,lb,ub,options)
% ipoptqp: solve QP with IpOpt
%
% [status,x,fopt] = ipoptqp_tanh(x0,tanh_spec,tanh_pairs_spec,Q,k,A,b,Aeq,beq,lb,ub,options)
%
% tanh_spec and tanh_pairs_spec must be matrices with 4 columns specifying
% [i s b w] for the expressions w*tanh(s*x[i]+b) (indexing starts with 1).
% In tanh_pairs_spec, two subsequent rows specify a pair
% w1*w2*tanh(s1*x[i1]+b1)*tanh(s2*x[i2]+b2).
%
% Options structure
%    options.maxabs           [double>0 (default Inf)]
%    options.boundstol        [double>0 (default 1e-8)]
%    options.extractbounds    [0 (off) or 1 (on) (default = 1, use to extract bounds from A, b, Aeq, Beq)]
%    options.qfactor          [double>=0 (default 2, use 2 for optimizing x'Qx instead of 0.5x'Qx)]
%    options.<ipopt-option>   set any of the IpOpt options
%
% please make sure to extract bounds if A or Aeq contain rows with a single
% nonzero element
%
%  status : 0 - optimal, 1 - infeasible, 2- max iter exceeded, 3- error

% written by Jan Poland, ABB.CHCRC.C1

if nargin<2 || size(tanh_spec,2)~=4
  help ipoptqp_tanh; return; 
end
if nargin<3, tanh_pairs_spec=[]; end
if ~isempty(tanh_pairs_spec) && (size(tanh_pairs_spec,2)~=4 || rem(size(tanh_pairs_spec,1),2)~=0)
  help ipoptqp_tanh; return; 
end
if nargin<4, Q = []; end
if nargin<5, k = []; end
if nargin<6, A = []; end
if nargin<7, b = []; end;
if nargin<8, Aeq = []; end;
if nargin<9, beq = []; end
if nargin<10, lb = []; end
if nargin<11, ub = []; end
if nargin<12, options = []; end

A = [Aeq;A];
b = [beq;b];
neq = length(beq);
nx = length(x0);
if isempty(Q), Q = sparse(nx,nx); end
if isempty(k), k = zeros(1,nx); end
  
extractbounds = 1; qfactor = 2; maxabs = inf; boundstol = 1e-8; printlevel = 0;
if isfield(options,'extractbounds')
  extractbounds = options.extractbounds; 
  options = rmfield(options,'extractbounds');
end
if isfield(options,'qfactor')
  qfactor = options.qfactor;
  options = rmfield(options,'qfactor');
end
if isfield(options,'maxabs')
  maxabs = options.maxabs;
  options = rmfield(options,'maxabs');
end
if isfield(options,'boundstol')
  boundstol = options.boundstol; 
  options = rmfield(options,'boundstol');
end
if isfield(options,'printlevel')
  printlevel = options.printlevel; 
  options = rmfield(options,'printlevel');
end
if isfield(options,'print_level')
  printlevel = options.print_level; 
  options = rmfield(options,'print_level');
end

if isempty(lb), lb=ones(1,nx)*(-maxabs); else lb=max(lb,-maxabs); end
if isempty(ub), ub=ones(1,nx)*(maxabs); else ub=min(ub,maxabs); end
lb = lb(:);
ub = ub(:);
x = zeros(size(lb));
jp = 1:length(x);
Q0 = []; k0 = [];

if extractbounds
  Q0 = Q; k0 = k;
  while 1
    nzA = sum(A~=0,2);
    ii = find(nzA==0);
    if ~isempty(ii)
      A(ii,:) = []; b(ii) = []; neq = neq-sum(ii<=neq);
    end
    
    ii = find(nzA==1);
    if ~isempty(ii)
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
      
      if any(lb>ub+boundstol)
        error('conflicting bounds!');
      end
      j = find(lb>ub);
      if ~isempty(j), lb(j)=ub(j); end
    end
    
    ii = find(lb>ub-boundstol);
    if isempty(ii), break; end
    
    x(jp(ii)) = (lb(ii)+ub(ii))/2;
    lb(ii) = []; ub(ii) = [];
    Q(ii,:) = []; Q(:,ii) = []; k(ii) = [];
    b = b - A(:,ii)*x(jp(ii));
    A(:,ii) = [];
    jp(ii) = [];
  end  
end

[iA,jA,vA]=find(A); iA=iA-1; jA=jA-1;

spQ = logical(Q);
tanh_spec = [tanh_pairs_spec;tanh_spec];
spQ(sub2ind(size(Q),tanh_spec(:,1),tanh_spec(:,1)))=true;
if ~isempty(tanh_pairs_spec)
  spQ(sub2ind(size(Q),tanh_pairs_spec(1:2:end,1),tanh_pairs_spec(2:2:end,1)))=true;
  spQ(sub2ind(size(Q),tanh_pairs_spec(2:2:end,1),tanh_pairs_spec(1:2:end,1)))=true;
end
[iQ,jQ]=find(tril(spQ)); 
vQ = full(Q(sub2ind(size(Q),iQ,jQ)));
JQ = sparse(size(Q)); JQ(sub2ind(size(Q),iQ,jQ)) = 1:length(iQ);
tanh_spec = [tanh_spec full(JQ(sub2ind(size(Q),tanh_spec(:,1),tanh_spec(:,1))))'-1];
if isempty(tanh_pairs_spec), jmixed=[]; else
  jmixed = repmat(max(full(JQ(sub2ind(size(Q),tanh_pairs_spec(1:2:end,1),tanh_pairs_spec(2:2:end,1)))),...
    full(JQ(sub2ind(size(Q),tanh_pairs_spec(2:2:end,1),tanh_pairs_spec(1:2:end,1))))),2,1);
end
tanh_spec = [tanh_spec [jmixed(:)-1;zeros(size(tanh_spec,1)-length(jmixed(:)),1)]];
tanh_spec(:,1)=tanh_spec(:,1)-1;  
iQ=iQ-1; jQ=jQ-1;

if qfactor~=1, vQ=qfactor*vQ; end

k = full(k(:))';
b = full(b(:))';
lb = full(lb(:)');
ub = full(ub(:)');

[status,x1,fopt] = ipoptqp_tanh_mex(lb,ub,k,iQ,jQ,vQ,b,iA,jA,vA,neq,...
  printlevel,x0,tanh_spec,size(tanh_pairs_spec,1),options);
x(jp) = x1;
if length(jp)<length(x)
  fopt = (qfactor/2)*x'*Q0*x+k0*x;
end

