function y = vmlSimplexProj(x)
% y = vmlSimplexProj(x)
% Project vector x onto the probability simplex.
% Algorithm from http://ttic.uchicago.edu/~wwang5/papers/SimplexProj.pdf
% JP / ABB.CH-RD.C1 / 2016-01-20

u = sort(x(:),1,'descend');
rho = find(u+(1-cumsum(u))./(1:length(u))'>0,1,'last');
lambda = 1/rho*(1-sum(u(1:rho)));
y = max(x+lambda,0);
