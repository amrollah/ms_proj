function [rmse] = critfun3(X,Y)
x=X(1:floor(3*size(Y,1)/4),:);
y=Y(1:floor(3*size(Y,1)/4),:);
x_t=X(floor(3*size(Y,1)/4)+1:end,:);
yt=Y(floor(3*size(Y,1)/4)+1:end,:);

svm_model = svmtrain(y, x, '-s 3 -t 2 -g 8 -c 250 -p 9');
[y_hat,Acc,~] = svmpredict(yt, x_t, svm_model);
rmse = sqrt(mean((yt - y_hat).^2));
disp(['Error svm regres: ' num2str(rmse)]);
end

