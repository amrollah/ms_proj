close all; 
clear all; 
% clc;
rng('default');
rng(2,'twister');
normal_regress=false;
svorm = false;
svm_regress = true;
remove_rare = false;
knn=false;

w=5;
fold=1;
img_save_path='';
prj_path='';
proj_path;
addpath(prj_path);
if ~exist('data','var')
    load('calc\clean_data_with_8cc_nan_corrected3.mat', 'data');
end
N = 30000; %9000; train size
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

% irr1 = cellfun(@(d) d.irr(1), data(train_ind));
% irr2 = cellfun(@(d) d.irr(2), data(train_ind));
% offset=45;
% clear_DNI = cellfun(@(d) d.clear_irr(2), data(train_ind));
% DNI = clear_DNI .* sun_flag;
% DHI1 = irr1 - cosd(zenith).*DNI;
% DHI2 = irr2 - max(0,cosd(zenith+offset)).*DNI;

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

% irr1_t = cellfun(@(d) d.irr(1), data(test_ind));
% irr2_t = cellfun(@(d) d.irr(2), data(test_ind));
% offset=45;
% clear_DNI_t = cellfun(@(d) d.clear_irr(2), data(test_ind));
% DNI_t = clear_DNI_t .* sun_flag_t;
% DHI1_t = irr1_t - cosd(zenith_t).*DNI_t;
% DHI2_t = irr2_t - max(0,cosd(zenith_t+offset)).*DNI_t;
% y=DHI2;
% yt=DHI2_t;
% sun_flag(sun_flag==0) = sat_fact(sun_flag==0);
x_img =[sun_flag;sat_fact;clouds;mean(cc_fact(1:4,:))];
x_no_img = [clear_diffuse;zenith];%;DAY;TM];
x_img=normalize(x_img,0,1);
x_no_img=normalize(x_no_img,0,1);
% sat_fact = max(sat_fact,sun_flag);

x = [
    sun_flag...
    ; sat_fact...
    ; cc_fact(5:8,:)...
    ; clear_diffuse...%.*(1-clouds./100)...
    ; zenith...
%     ; clouds...
    ; mean(cc_fact(1:4,:))...
%     ; DAY...
%     ; TM
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
x_no_img_t = [clear_diffuse_t;zenith_t];%;DAY_t;TM_t];
x_img_t=normalize(x_img_t,0,1);
x_no_img_t=normalize(x_no_img_t,0,1);
% sat_fact_t = max(sat_fact_t,sun_flag_t);

test = [
    sun_flag_t...
    ; sat_fact_t...
    ; cc_fact_t(5:8,:)...
    ; clear_diffuse_t...%.*(1-clouds_t./100)...
    ; zenith_t...
%     ; clouds_t...
    ; mean(cc_fact_t(1:4,:))...
%     ; DAY_t...
%     ; TM_t
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
y_hat = max(0,predict(model,x_t'));
rmse1 = sqrt(mean((yt - y_hat').^2));
disp(['Error linear regres: ' num2str(rmse1)]);
% x=x_no_img;
maxdev = chi2inv(.60,1);     
opt = statset('display','iter',...
              'TolFun',maxdev,...
              'TolTypeFun','abs');

inmodel = sequentialfs(@critfun2,[x',sqrt(x'),(x').^2],y',...
                       'cv','none',...
                       'nullmodel',false,...
                       'options',opt,...
                       'direction','forward');


new_x = [x',sqrt(x'),(x').^2];
model = fitlm(new_x(:,inmodel),y','interactions');
% x_t=x_no_img_t;
new_x_t = [x_t',sqrt(x_t'),(x_t').^2];
y_hat = max(0,predict(model,new_x_t(:,inmodel)));
new_rmse1 = sqrt(mean((yt - y_hat').^2));
disp(['Error linear regres: ' num2str(new_rmse1)]);
err1 = yt - y_hat';

y_hat_tr = max(0,predict(model,new_x(:,inmodel)));
disp(['Error linear regres: ' num2str(sqrt(mean((y - y_hat_tr').^2)))]);
err3 = y - y_hat_tr';

maxdev = chi2inv(.60,1);     
opt = statset('display','iter',...
              'TolFun',maxdev,...
              'TolTypeFun','abs');

inmodel = sequentialfs(@critfun2,[x_no_img',sqrt(x_no_img'),(x_no_img').^2],y',...
                       'cv','none',...
                       'nullmodel',false,...
                       'options',opt,...
                       'direction','forward');


new_x = [x_no_img',sqrt(x_no_img'),(x_no_img').^2];
model = fitlm(new_x(:,inmodel),y','interactions');
% x_t=x_no_img_t;
new_x_t = [x_no_img_t',sqrt(x_no_img_t'),(x_no_img_t').^2];
y_hat = max(0,predict(model,new_x_t(:,inmodel)));
new_rmse1 = sqrt(mean((yt - y_hat').^2));
disp(['Error linear regres: ' num2str(new_rmse1)]);
err2 = yt - y_hat';

y_hat_tr = max(0,predict(model,new_x(:,inmodel)));
disp(['Error linear regres: ' num2str(sqrt(mean((y - y_hat_tr').^2)))]);
err4 = y - y_hat_tr';



figure(101);
scatter(y_hat,yt,'b','.');
lsline;
hold on;
grid on;
plot([0 400], [0 400],'r-');
xlabel('Predit Irradiance W/M^2');
ylabel('Measured Irradiance W/M^2');
xlim([0 400]);
ylim([0 400]);
title('Standard regression (non-image feature)');
    
end

% [ranked,weights] = relieff(x',y',3);
% figure(919);
% bar(weights(ranked));
% xlabel('Predictor rank');
% ylabel('Predictor importance weight');
if knn
%% k-nearest neighbors
maxdev = chi2inv(.60,1);     
opt = statset('display','iter',...
              'TolFun',maxdev,...
              'TolTypeFun','abs');

inmodel = sequentialfs(@critfun,x',y',...
                       'cv','none',...
                       'nullmodel',false,...
                       'options',opt,...
                       'direction','forward');
                   
% [IDX,D] = knnsearch(x([1 6 8 9 11 12],:)',x_t([1 6 8 9 11 12],:)','K',2);
[IDX,D] = knnsearch(x',x_t','K',10,'Distance','minkowski','P',1);
WD=D./repmat(sum(D,2),[1,10]);
knn_y_hat = sum(y(IDX).*WD,2);
% err5 = abs(yt - knn_y_hat').*(log(yt)/log(100));
% disp(['Error knn: ' num2str(mean(err5)), '   std: ', num2str(std(err5))]);
rmse5 = sqrt(mean((yt - knn_y_hat').^2));
% rmse = sqrt((yt - knn_y_hat').^2);
% knn_R_sq = 1-sum((yt - knn_y_hat').^2)/sum(yt.^2);
% knn_MBE = mean(yt - knn_y_hat');
disp(['Error knn regres: ' num2str(rmse5)]);
% result_show(data(test_ind),knn_y_hat',yt,IDX,data(train_ind));
err1=yt - knn_y_hat';

[IDX,D] = knnsearch(x',x','K',10,'Distance','minkowski','P',1);
WD=D./repmat(sum(D,2),[1,10]);
knn_y_hat = sum(y(IDX).*WD,2);
rmse5 = sqrt(mean((y - knn_y_hat').^2));
disp(['Error knn regres: ' num2str(rmse5)]);
err3=y - knn_y_hat';

[IDX,D] = knnsearch(x_no_img',x_no_img_t','K',10,'Distance','minkowski','P',1);
WD=D./repmat(sum(D,2),[1,10]);
knn_y_hat = sum(y(IDX).*WD,2);
rmse = sqrt(mean((yt - knn_y_hat').^2));
disp(['Error knn regres: ' num2str(rmse)]);
err2=yt - knn_y_hat';

[IDX,D] = knnsearch(x_no_img',x_no_img','K',10,'Distance','minkowski','P',1);
WD=D./repmat(sum(D,2),[1,10]);
knn_y_hat = sum(y(IDX).*WD,2);
rmse = sqrt(mean((y - knn_y_hat').^2));
disp(['Error knn regres: ' num2str(rmse)]);
err4=y - knn_y_hat';

[y1,x1] = hist(err1,50); [y2,x2] = hist(err2,50); [y3,x3] = hist(err3,50); [y4,x4] = hist(err4,50);
figure(13); h1=plot(x1,y1,'b'); grid on; hold on; h2=plot(x2,y2,'r'); h3=plot(x3,y3,'b--'); h4=plot(x4,y4,'r--'); xlim([-300,300]); legend([h1,h2,h3,h4],{'All features(test)'; 'non-image features(test)';'All features(train)';'non-image features(train)'});
xlabel('Error (W/m^2)'); ylabel('log(frequency)'); title('Histogram of errors'); hold off;


figure(2);
% values = hist3([knn_y_hat(:)./clear_diffuse_t' yt(:)./clear_diffuse_t'],[20 20]);
% newmap = jet;
% newmap(1,:) = [1 1 1];
% colormap(newmap); 
plot(knn_y_hat, yt,'.');
% imagesc(values)
% lowrange = 1;
% highrange = 10;
% caxis manual
% caxis([lowrange highrange]);
% colorbar
% axis equal
% axis xy
ylabel('DHI observed (W/m^2)');
xlabel('DHI estimated (W/m^2)');
hold on;
plot([0 400], [0 400],'r-');
title('Correlation of DHI estimation from K-NN regresor');

% stairs
% hist
[y1,x1] = hist(err1,50); [y2,x2] = hist(err2,50); [y3,x3] = hist(err3,50); [y4,x4] = hist(err4,50);
figure(10); h1=stairs(x1,y1,'b'); grid on; hold on; h2=stairs(x2,y2,'r'); h3=stairs(x3,y3,'b--'); h4=stairs(x4,y4,'r--'); xlim([-300,300]); legend([h1,h2,h3,h4],{'All features(test)'; 'non-image features(test)';'All features(train)';'non-image features(train)'});
xlabel('Error (W/m^2)'); ylabel('log(frequency)'); title('Histogram of errors'); hold off;

[y1,x1] = hist(err1,50); [y2,x2] = hist(err2,50); [y3,x3] = hist(err3,50); [y4,x4] = hist(err4,50);
figure(12); h1=plot(x1,y1,'b'); grid on; hold on; h2=plot(x2,y2,'r'); h3=plot(x3,y3,'b--'); h4=plot(x4,y4,'r--'); xlim([-300,300]); legend([h1,h2,h3,h4],{'All features(test)'; 'non-image features(test)';'All features(train)';'non-image features(train)'});
xlabel('Error (W/m^2)'); ylabel('log(frequency)'); title('Histogram of errors'); hold off;

end                   
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
% dlmwrite('E:/ABB/svorim/d_train.0',train,' ');
% dlmwrite('E:/ABB/svorim/d_test.0',[test;yt]',' ');


% Another library for svm-regression
% svrobj = svr_trainer(x',y',400,0.000000025,'gaussian',0.5);
% y_hat = svrobj.predict(test);


if svm_regress
% libSVM regressor
maxdev = chi2inv(.60,1);     
opt = statset('display','iter',...
              'TolFun',maxdev,...
              'TolTypeFun','abs');

inmodel = sequentialfs(@critfun3,x',y',...
                       'cv','none',...
                       'nullmodel',false,...
                       'options',opt,...
                       'direction','forward');
                   
svm_model = svmtrain(y', x', '-s 3 -t 2 -g 8 -c 250 -p 9');
[y_hat,Acc,~] = svmpredict(yt', test', svm_model);
rmse4 = sqrt(mean((yt - y_hat').^2));
disp(['Error svm regres: ' num2str(rmse4)]);
err1 = yt - y_hat';

[y_hat,Acc,~] = svmpredict(y', x', svm_model);
rmse4 = sqrt(mean((y - y_hat').^2));
disp(['Error svm regres: ' num2str(rmse4)]);
err3 = y - y_hat';

svm_model = svmtrain(y', x_no_img', '-s 3 -t 2 -g 8 -c 250 -p 9');
[y_hat,Acc,~] = svmpredict(yt', x_no_img_t', svm_model);
rmse4 = sqrt(mean((yt - y_hat').^2));
disp(['Error svm regres: ' num2str(rmse4)]);
err2 = yt - y_hat';

[y_hat,Acc,~] = svmpredict(y', x_no_img', svm_model);
rmse4 = sqrt(mean((y - y_hat').^2));
disp(['Error svm regres: ' num2str(rmse4)]);
err4 = y - y_hat';


[y1,x1] = hist(err1,50); [y2,x2] = hist(err2,50); [y3,x3] = hist(err3,50); [y4,x4] = hist(err4,50);
figure(14); h1=plot(x1,y1,'b'); grid on; hold on; h2=plot(x2,y2,'r'); h3=plot(x3,y3,'b--'); h4=plot(x4,y4,'r--'); xlim([-300,300]); legend([h1,h2,h3,h4],{'All features(test)'; 'non-image features(test)';'All features(train)';'non-image features(train)'});
xlabel('Error (W/m^2)'); ylabel('frequency'); title('Histogram of errors'); hold off;

figure(108);
scatter(y_hat,yt,'b','.');
lsline;
hold on;
grid on;
plot([0 400], [0 400],'r-');
xlabel('Predit Irradiance W/M^2');
ylabel('Measured Irradiance W/M^2');
xlim([0 400]);
ylim([0 400]);
title('SVR result (all the features)');

result_show(data(test_ind),y_hat',yt,IDX,data(train_ind));

err4 = abs(yt - y_hat').*(log(yt)/log(100));
err14 = abs(yt - y_hat')./log(yt);

R_sq = 1-sum((yt - y_hat').^2)/sum(yt.^2);
MBE = mean(yt - y_hat');
disp(['Error svm regres: ' num2str(rmse4), '   std: ', num2str(std(err4))]);
end