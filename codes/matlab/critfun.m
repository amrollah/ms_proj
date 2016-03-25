function [rmse] = critfun(X,Y)
x=X(1:floor(3*size(Y,1)/4),:);
y=Y(1:floor(3*size(Y,1)/4),:);
x_t=X(floor(3*size(Y,1)/4)+1:end,:);
yt=Y(floor(3*size(Y,1)/4)+1:end,:);
K=3;
[IDX,D] = knnsearch(x,x_t,'K',K);
sm = sum(D,2);
WD=D./repmat(sm,[1,K]);
WD(isnan(WD)==1)=1/K;
knn_y_hat = sum(y(IDX).*WD,2);
rmse = sqrt(mean((yt - knn_y_hat).^2));
disp(['Error knn regres: ' num2str(rmse)]);
end

