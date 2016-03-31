function [rmse] = critfun2(X,Y)
x=X(1:floor(3*size(Y,1)/4),:);
y=Y(1:floor(3*size(Y,1)/4),:);
x_t=X(floor(3*size(Y,1)/4)+1:end,:);
yt=Y(floor(3*size(Y,1)/4)+1:end,:);
model = fitlm(x,y,'interactions');
y_hat = max(0,predict(model,x_t));
rmse = sqrt(mean((yt - y_hat).^2));
disp(['Error knn regres: ' num2str(rmse)]);
end

