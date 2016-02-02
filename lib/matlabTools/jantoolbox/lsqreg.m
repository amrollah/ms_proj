function b = lsqreg(X,y,lambda)
if all(size(lambda)==1), lambda = lambda*eye(size(X,2));
elseif diff(size(lambda))~=0, lambda=diag(lambda); 
end 
[Q,R]=qr(X'*X+lambda);
b = R\(Q'*(X'*y));
