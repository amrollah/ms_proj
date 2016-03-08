close all; 
% clear all; 
clc;
w=10;
fold=2;
img_save_path='';
prj_path='';
proj_path;
addpath(prj_path);
if ~exist('data','var')
    load('calc\clean_data_with_8cc_nan_corrected.mat', 'data');
end
A=1:length(data);
train_ind = 1:50:length(data);
test_ind = A(~ismember(A,train_ind));
y = cellfun(@(d) d.corr_tilt_diff, data(train_ind));
% round_y = int64(y);
% labels = unique(round_y);
% target = arrayfun(@(yi) find(labels==yi), round_y);
% y=uint8(normalize(y,1,100));

round_y = w*ceil(y/w);
labels = unique(round_y);
target = arrayfun(@(yi) find(labels==yi), round_y);
nlabels=labels;
for i=1:length(labels)
    trg = find(target==i);
    if length(trg)<fold
        nlabels(nlabels==labels(i)) = [];
        train_ind(trg)=[];
    end
end
if length(labels)>length(nlabels)
    y = cellfun(@(d) d.corr_tilt_diff, data(train_ind));
    round_y = w*ceil(y/w);
    labels = unique(round_y);
    target = arrayfun(@(yi) find(labels==yi), round_y);
end
sun_flag = cellfun(@(d) d.sun_flag, data(train_ind));
azimuth = cellfun(@(d) d.azimuth, data(train_ind));
zenith = cellfun(@(d) d.zenith, data(train_ind));
clouds  = cellfun(@(d) d.clouds, data(train_ind));
sat_fact = cellfun(@(d) d.sat_fact, data(train_ind));
cc_fact = cell2mat(cellfun(@(d) d.cc_fact', data(train_ind),'UniformOutput', false));
clear_diffuse = cellfun(@(d) d.clear_irr(3), data(train_ind));
time = cellfun(@(d) d.time, data(train_ind));
DV  = datevec(time);
TM = (DV(:,4)*60+DV(:,5))'; 
DV  = DV(:, 1:3); DV2 = DV;
DV2(:, 2:3) = repmat([7,1],[size(DV2,1),1]);
DAY = double(abs(int64(time' - datenum(DV2))))';

x = [sun_flag; sat_fact; clear_diffuse; zenith; cc_fact];
x=normalize(x,0,1);
train = [x; target]';

yt = cellfun(@(d) d.corr_tilt_diff, data(test_ind));
sun_flag_t = cellfun(@(d) d.sun_flag, data(test_ind));
azimuth_t = cellfun(@(d) d.azimuth, data(test_ind));
zenith_t = cellfun(@(d) d.zenith, data(test_ind));
clouds_t  = cellfun(@(d) d.clouds, data(test_ind));
cc_fact_t = cell2mat(cellfun(@(d) d.cc_fact', data(test_ind),'UniformOutput', false));
sat_fact_t = cellfun(@(d) d.sat_fact, data(test_ind));
clear_diffuse_t = cellfun(@(d) d.clear_irr(3), data(test_ind));
time_t = cellfun(@(d) d.time, data(test_ind));
DV  = datevec(time_t);
TM_t = (DV(:,4)*60+DV(:,5))'; 
DV  = DV(:, 1:3); DV2 = DV;
DV2(:, 2:3) = repmat([7,1],[size(DV2,1),1]);
DAY_t = double(abs(int64(time_t' - datenum(DV2))))';

test = [sun_flag_t; sat_fact_t; clear_diffuse_t; zenith_t; cc_fact_t]';
test=normalize(test,0,1);

model = fitlm(x',y','interactions');
y_hat = max(0,predict(model,test));
figure;
scatter(y_hat,yt,'b','.');
lsline;
xlabel('predit');
ylabel('measured');
grid on;
title('Standard regression');


% dlmwrite('E:/ABB/svorim/d_train.0',train,' ');
% dlmwrite('E:/ABB/svorim/d_test.0',test,' ');
% 

% Another library for svm-regression
% svrobj = svr_trainer(x',y',400,0.000000025,'gaussian',0.5);
% y_hat = svrobj.predict(test);

% yet another SVM regressor
% svm_model = svmtrain(y', x', '-s 3 -t 0 -c 20 -p 1');
% [y_hat,Acc,~] = svmpredict(yt', test, svm_model);

m = svorim(x,target);
target_t= m(test');
figure;
scatter(target_t,yt,'b','.');
g_ln = lsline;
xlabel('predit');
ylabel('measured');
title('SV ordinal regression machine');
grid on;
hold on;

grid_sz = 30;
dxi = ceil(length(test)/grid_sz);
xi = 1:dxi:length(test);

xd = floor(length(test)/grid_sz);
ln = {};
for k=2:grid_sz
    xk = target_t(xi(k-1):xi(k));
    yk = yt(xi(k-1):xi(k));
    p = polyfit(xk,yk,1);
    plot(xk,polyval(p,xk));
    ln{k-1} = p;
    hold on;
end



dxi = ceil(length(target_t)/30);
nodes = 1:dxi:length(target_t);
nnodes = length(nodes);
scale = 3;
xi = target_t(nodes);
dx=mean(xi(2:end) - xi(1:end-1));
dm = scale *dx * ones(1, nnodes);
[PHI, DPHI, DDPHI] = MLS1DShape(1, nnodes, xi, length(target_t), target_t, dm, 'GAUSS', .3);
% Curve fitting
yi  = yt(nodes);    % Nodal function values
yh  = PHI * yi';  % Approximation function
err = norm(yt' - yh) / norm(yt) * 100;  % Relative error norm in approximation function


figure;
scatter(yh',yt,'b','.');
lsline;
xlabel('predit');
ylabel('measured');
grid on;
title('SVM regression');
