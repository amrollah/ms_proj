close all; clc;
w=10;
fold=2;
img_save_path='';
prj_path='';
proj_path;
addpath(prj_path);
if ~exist('data','var')
    load('calc\data_with_sat_fact.mat', 'data');
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
% nlabels=labels;
% for i=1:length(labels)
%     trg = find(target==i);
%     if length(trg)<fold
%         nlabels(nlabels==labels(i)) = [];
%         train_ind(trg)=[];
%     end
% end
% if length(labels)>length(nlabels)
%     y = cellfun(@(d) d.corr_tilt_diff, data(train_ind));
%     round_y = w*ceil(y/w);
%     labels = unique(round_y);
%     target = arrayfun(@(yi) find(labels==yi), round_y);
% end
sun_flag = cellfun(@(d) d.sun_flag, data(train_ind));
azimuth = cellfun(@(d) d.azimuth, data(train_ind));
zenith = cellfun(@(d) d.zenith, data(train_ind));
clouds  = cellfun(@(d) d.clouds, data(train_ind));
sat_fact = cellfun(@(d) d.sat_fact, data(train_ind));
clear_diffuse = cellfun(@(d) d.clear_irr(3), data(train_ind));
time = cellfun(@(d) d.time, data(train_ind));
DV  = datevec(time);
TM = (DV(:,4)*60+DV(:,5))'; 
DV  = DV(:, 1:3); DV2 = DV;
DV2(:, 2:3) = repmat([7,1],[size(DV2,1),1]);
DAY = double(abs(int64(time' - datenum(DV2))))';

x = [sun_flag; sat_fact; clear_diffuse; zenith; clouds];
x=normalize(x,0,1);
train = [sun_flag; sat_fact; clear_diffuse; zenith; clouds; target]';

yt = cellfun(@(d) int64(d.corr_tilt_diff), data(test_ind));
sun_flag_t = cellfun(@(d) d.sun_flag, data(test_ind));
azimuth_t = cellfun(@(d) d.azimuth, data(test_ind));
zenith_t = cellfun(@(d) d.zenith, data(test_ind));
clouds_t  = cellfun(@(d) d.clouds, data(test_ind));
sat_fact_t = cellfun(@(d) d.sat_fact, data(test_ind));
clear_diffuse_t = cellfun(@(d) d.clear_irr(3), data(test_ind));
time_t = cellfun(@(d) d.time, data(test_ind));
DV  = datevec(time_t);
TM_t = (DV(:,4)*60+DV(:,5))'; 
DV  = DV(:, 1:3); DV2 = DV;
DV2(:, 2:3) = repmat([7,1],[size(DV2,1),1]);
DAY_t = double(abs(int64(time_t' - datenum(DV2))))';

test = [sun_flag_t; sat_fact_t; clear_diffuse_t; zenith_t; clouds_t]';
test=normalize(test,0,1);

model = fitlm(x',y','quadratic');

disp_ind = 1:500:length(test);
est_diff = predict(model,test);
figure(210);
h(1) = plot(time_t(disp_ind),est_diff(disp_ind),'r.-');
hold on;
h(2) = plot(time_t(disp_ind),yt(disp_ind),'b.-');
legend([h(1),h(2)],{'Predict'; 'Observed'});
ylabel('Diffuse');
grid on;
datetickzoom;
anova(model)

% dlmwrite('E:/ABB/svorim/d_train.0',train,' ');
% dlmwrite('E:/ABB/svorim/d_test.0',test,' ');

m = svorim(x,target);
target_t= m(test');
figure(211);
h(1) = plot(time_t,labels(ceil(target_t)),'r.-');
hold on;
h(2) = plot(time_t,yt,'b.-');
legend([h(1),h(2)],{'Predict'; 'Observed'});
ylabel('Diffuse');
grid on;
datetickzoom;