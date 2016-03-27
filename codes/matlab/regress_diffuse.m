close all; 
clear all; 
% clc;
rng('default');
rng(2,'twister');
normal_regress=false;
svorm = false;
svm_regress = true;
remove_rare = false;
w=5;
fold=1;
img_save_path='';
prj_path='';
proj_path;
addpath(prj_path);
if ~exist('data','var')
    load('calc\clean_data_with_8cc_nan_corrected3.mat', 'data');
end
N = 30000; % train size
train_ind = [];
A=1:length(data);
all_y = cellfun(@(d) d.corr_tilt_diff, data);
bin_100 = find(all_y<100);
train_ind = [train_ind bin_100(randperm(length(bin_100),floor(N/5)))];
bin_200 = find(all_y>=100 & all_y<200);
train_ind = [train_ind bin_200(randperm(length(bin_200),floor(N/5)))];
bin_300 = find(all_y>=200 & all_y<300);
train_ind = [train_ind bin_300(randperm(length(bin_300),floor(min(N/5,.7*length(bin_300)))))];
bin_400 = find(all_y>=300 & all_y<400);
train_ind = [train_ind bin_400(randperm(length(bin_400),floor(min(N/5,.7*length(bin_400)))))];
bin_500 = find(all_y>=400 & all_y<600);
train_ind = [train_ind bin_500(randperm(length(bin_500),floor(min(N/10,.7*length(bin_500)))))];

% train_ind = randperm(length(data),500);
% while true
test_ind = (~ismember(A,train_ind));
y = cellfun(@(d) d.corr_tilt_diff, data(train_ind));
round_y = w*ceil(y/w);
labels = unique(round_y);
target = arrayfun(@(yi) find(labels==yi), round_y);
% removing the very rare data, probabely outliers
if remove_rare
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
end
%% extracting training and test data 
sun_flag = cellfun(@(d) d.sun_flag, data(train_ind));
sun_flag=normalize(sun_flag,0,1);
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

irr1 = cellfun(@(d) d.irr(1), data(train_ind));
irr2 = cellfun(@(d) d.irr(2), data(train_ind));
offset=45;
clear_DNI = cellfun(@(d) d.clear_irr(2), data(train_ind));
DNI = clear_DNI .* sun_flag;
DHI1 = irr1 - cosd(zenith).*DNI;
DHI2 = irr2 - max(0,cosd(zenith+offset)).*DNI;

%test features
yt = cellfun(@(d) d.corr_tilt_diff, data(test_ind));
sun_flag_t = cellfun(@(d) d.sun_flag, data(test_ind));
sun_flag_t=normalize(sun_flag_t,0,1);
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

irr1_t = cellfun(@(d) d.irr(1), data(test_ind));
irr2_t = cellfun(@(d) d.irr(2), data(test_ind));
offset=45;
clear_DNI_t = cellfun(@(d) d.clear_irr(2), data(test_ind));
DNI_t = clear_DNI_t .* sun_flag_t;
DHI1_t = irr1_t - cosd(zenith_t).*DNI_t;
DHI2_t = irr2_t - max(0,cosd(zenith_t+offset)).*DNI_t;
y=DHI2;
yt=DHI2_t;

x_img =[sun_flag;sat_fact;clouds;mean(cc_fact(1:4,:))];
x_no_img = [clear_diffuse;zenith;DAY;TM];
x_img=normalize(x_img,0,1);
x_no_img=normalize(x_no_img,0,1);

x = [
    sun_flag...
    ; sat_fact...
    ; clear_diffuse...%.*(1-clouds./100)...
    ; zenith...
    ; clouds...
    ; mean(cc_fact(1:4,:))...
    ; DAY...
    ; TM
    ];
x = [
    x...
%     ;x.^2 ...
%     ;x.^.5 ...
%     ;sat_fact.*(1./(sun_flag+.1))...
%     ;(1./(sun_flag+.1)).*clouds;...
%     ;sat_fact.*clear_diffuse;sat_fact.*clouds;...
%     ;clear_diffuse./max(0.1,clouds);...
%     ;zenith.*clouds;zenith.*DAY;zenith.*TM...
%     ;clouds.*DAY;clouds.*TM...
%     ;DAY.*TM...
    ];
x=normalize(x,0,1);
train = [x; y]';

x_img_t =[sun_flag_t;sat_fact_t;clouds_t;mean(cc_fact_t(1:4,:))];
x_no_img_t = [clear_diffuse_t;zenith_t;DAY_t;TM_t];
x_img_t=normalize(x_img_t,0,1);
x_no_img_t=normalize(x_no_img_t,0,1);

test = [
    sun_flag_t...
    ; sat_fact_t...
    ; clear_diffuse_t...%.*(1-clouds_t./100)...
    ; zenith_t...
    ; clouds_t...
    ; mean(cc_fact_t(1:4,:))...
    ; DAY_t...
    ; TM_t
    ];
test = [
    test...
%      ;test.^2 ...
%      ;test.^.5 ...
%     ;sat_fact_t.*(1./(sun_flag_t+.1))...
%     ;(1./(sun_flag_t+.1)).*clouds_t;...
%     ;sat_fact_t.*clear_diffuse_t;sat_fact_t.*clouds_t;...
%     ;clear_diffuse_t./max(0.1,clouds_t);...
%     ;zenith_t.*clouds_t;zenith_t.*DAY_t;zenith_t.*TM_t...
%     ;clouds_t.*DAY_t;clouds_t.*TM_t...
%     ;DAY_t.*TM_t...
    ];
test=normalize(test,0,1);
x_t = test;

%% Standard regression
if normal_regress
model = fitlm(x',y','interactions');
y_hat = max(0,predict(model,test'));
figure(100);
scatter(y_hat,yt,'b','.');
lsline;
xlabel('predit');
ylabel('measured');
grid on;
title('Standard regression');
end

[ranked,weights] = relieff(x',y',3);
figure(919);
bar(weights(ranked));
xlabel('Predictor rank');
ylabel('Predictor importance weight');

%% k-nearest neighbors
[IDX,D] = knnsearch(x([1,3,4,5,6,7,8],:)',x_t([1,3,4,5,6,7,8],:)','K',2);
WD=D./repmat(sum(D,2),[1,2]);
knn_y_hat = sum(y(IDX).*WD,2);
err5 = abs(yt - knn_y_hat').*(log(yt)/log(100));
% disp(['Error knn: ' num2str(mean(err5)), '   std: ', num2str(std(err5))]);
rmse5 = sqrt(mean((yt - knn_y_hat').^2));
rmse = sqrt((yt - knn_y_hat').^2);
knn_R_sq = 1-sum((yt - knn_y_hat').^2)/sum(yt.^2);
knn_MBE = mean(yt - knn_y_hat');
disp(['Error knn regres: ' num2str(rmse5)]);
result_show(data(test_ind),knn_y_hat',yt,IDX,data(train_ind));

maxdev = chi2inv(.70,1);     
opt = statset('display','iter',...
              'TolFun',maxdev,...
              'TolTypeFun','abs');

inmodel = sequentialfs(@critfun,x',y',...
                       'cv','none',...
                       'nullmodel',false,...
                       'options',opt,...
                       'direction','forward');
                   
% cross-validated knn model
% mdl = fitcknn(x',y ,'NumNeighbors',2);
% knn_y_hat = predict(mdl,x_t');
% rmse5 = sqrt(mean((yt - knn_y_hat').^2));

%% SV Ordinal regression
if svorm
kernel = 2;
c=1000;
range=[labels(1),labels(end)];
m = svorim(x,target,kernel,c,range);
target_t= m(test);
% figure(1);
% scatter(target_t,yt,'b','.');
% xlim([0 400]);
% ylim([0 400]);
% lsline;
% hold on;
g_ln = polyfit(target_t,yt,1);
% yh=polyval(g_ln,target_t);
% plot(target_t,yh,'c.','LineWidth',1);
% xlabel('predit');
% ylabel('measured');
% title('SV ordinal regression machine');
% grid on;
% hold on;

%% smoothing
xi = floor(min(target_t)):20:ceil(max(target_t));
grid_sz = length(xi);
ln = {};
y_reg = [];
x_reg = [];
for k=2:grid_sz
    xki = find(target_t>=xi(k-1) & target_t<xi(k));
    xkv = target_t(xki);
    if length(xki)<300
        xkv=sort(xkv);
        yk=polyval(g_ln,xkv);
        ln{end+1} = g_ln;
    else 
        yk = yt(xki);
        p = polyfit(xkv,yk,1);
        xkv=sort(xkv);
        if p(1)>g_ln(1)
            yk=polyval(g_ln,xkv);
            ln{end+1} = g_ln;
        else
            yk=polyval(p,xkv);
            ln{end+1} = p;
        end
    end
%     plot(xkv,yk);
    x_reg = [x_reg,xkv];
    y_reg = [y_reg,yk];
    hold on;
end
% y_reg = medfilt1(y_reg);
y_smooth = smooth(x_reg,y_reg,0.05,'moving')';
% plot(x_reg,y_smooth,'r-','LineWidth',4);
% save('calc\final_reg_model.mat','x_reg','y_smooth','ln');

%% QUERY
yhat = interp1(x_reg,y_smooth,target_t,'spline');
% figure(2);plot(yhat,yt,'b.','LineWidth',3);
% xlabel('Predicted');
% ylabel('Measured');
% xlim([0 400]);
% ylim([0 400]);
% hold on;
% plot([0 400], [0 400],'r-');
 
% err1 = abs(yt - yhat).*(yt./100);  % Error
err1 = abs(yt - yhat).*(log(yt)/log(100)); 
err2 = abs(yt - target_t).*(log(yt)/log(100));
err3 = abs(yt - y_hat').*(log(yt)/log(100));
% RMSE2 = sqrt(mean((yt - target_t).^2));
disp(['RMSE interp: ' num2str(mean(err1)), '   std: ', num2str(std(err1))]);
disp(['RMSE svorim: ' num2str(mean(err2)), '   std: ', num2str(std(err2))]);
disp(['RMSE regres: ' num2str(mean(err3)), '   std: ', num2str(std(err3))]);

result_show(data(test_ind),target_t,yt);
end
% large_errs = find(err2>40);
% trn_ind = test_ind(large_errs(randperm(length(large_errs),400)));
% train_ind = [train_ind,trn_ind];
% pause();
% end

%% Other Methods
dlmwrite('E:/ABB/svorim/d_train.0',train,' ');
dlmwrite('E:/ABB/svorim/d_test.0',[test;yt]',' ');


% Another library for svm-regression
% svrobj = svr_trainer(x',y',400,0.000000025,'gaussian',0.5);
% y_hat = svrobj.predict(test);


if svm_regress
% libSVM regressor
svm_model = svmtrain(y', x(:,:)', '-s 3 -t 2 -g 8 -c 250 -p 9');
% cross_valid
[y_hat,Acc,~] = svmpredict(yt', test(:,:)', svm_model);

result_show(data(test_ind),y_hat',yt);

err4 = abs(yt - y_hat').*(log(yt)/log(100));
err14 = abs(yt - y_hat')./log(yt);
rmse4 = sqrt(mean((yt - y_hat').^2));
R_sq = 1-sum((yt - y_hat').^2)/sum(yt.^2);
MBE = mean(yt - y_hat');
disp(['Error svm regres: ' num2str(rmse4), '   std: ', num2str(std(err4))]);
end