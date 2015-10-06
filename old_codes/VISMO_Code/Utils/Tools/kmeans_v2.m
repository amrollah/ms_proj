% Simplified fast k-means implementation. 
% Burak Zeydan 21. Aug. 2014

function [idx, C, i] = kmeans_v2( X, k )
maxiter = 100;
prs = randi([1,size(X,2)],1,k);
C_s = zeros(size(X));
C_s(:,prs) = 1;
C = X(logical(C_s));
idx = zeros(1,size(X,2));
i = 0;
for i = 1:maxiter
    idx_pr = idx;
    X_cmp = (repmat(X,k,1)-repmat(C(:),1,size(X,2))).^2;
    [~,idx] = min(reshape(sum(reshape(X_cmp,size(X_cmp,1)/k,[],k),1),k,[]),[],1);
    for j = 1:k
        C(:,j) = mean(X(:,idx==j),2);
    end
    if(sum(idx == idx_pr)>0.98*size(X,2))
        break;
    end
end

end

