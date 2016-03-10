
%Usage: ./fselect.py training_file [testing_file]
%Output files: .fscore shows importance of features, .select gives the running log, and .pred gives testing results.

%If you could export your data in libsvm format I would suggest you to run easy.py tool on it and estimate the parameters. Then you could set the parameters (or even import the created model) in your matlab-libsvm environment.

%I could achieve accuracy improvements from <70% to >90% without understanding much of tuning SVM only by using this tool and reading this doc.

%# grid of parameters
folds = 5;
[eps,gamma] = meshgrid(0:5, -10:2:7);
C = max(y) - min(y);

%# grid search, and cross-validation
cv_acc = zeros(numel(gamma),1);
for i=1:numel(gamma)
    cv_acc(i) = svmtrain(y', x', ...
                    sprintf('-c %f -g %f -p %f -v %d', C, 2^gamma(i), eps(i), folds));
end

%# pair (C,gamma) with best accuracy
[~,idx] = max(cv_acc);

%# contour plot of paramter selection
contour(eps, gamma, reshape(cv_acc,size(eps))), colorbar
hold on
plot(eps(idx), gamma(idx), 'rx')
text(eps(idx), gamma(idx), sprintf('Acc = %.2f %%',cv_acc(idx)), ...
    'HorizontalAlign','left', 'VerticalAlign','top')
hold off
xlabel('log_2(eps)'), ylabel('log_2(\gamma)'), title('Cross-Validation Accuracy')

%# now you can train you model using best_C and best_gamma
best_eps = eps(idx);
best_gamma = 2^gamma(idx);

svm_model = svmtrain(y', x', sprintf('-c %f -g %f -p %f', C, best_gamma, best_eps));
					