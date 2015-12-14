% close all; 
clear all; 
% clc;

%% settings
addpath(genpath('C:\Users\chamsei\Documents\GitHub\ms_proj\old_codes\vismolib'))
addpath(genpath('C:\Users\chamsei\Documents\GitHub\ms_proj\old_codes\solvers'));
vmlconfig_cavriglia;
conf = evalin('base','VMLCONF');

days = {'2015_07_19',
        '2015_08_03',
        '2015_10_24',
        '2015_11_12'};

test_days={'2015_11_12',};
load_model = false;

if load_model
    load('lnModel.mat');
else
all_Irr=[];
all_temp=[];
all_sun_zth=[];
all_pw=[];
for d=1:numel(days)
    date = char(days(d));
    obj = vmlSeq(date);

    Irr_ts = timeseries(obj.Irr(:,2), obj.Irr(:,1));
    Irr = resample(Irr_ts, obj.P(:,1)); 
    all_Irr = [all_Irr;Irr.data];
    Temp_ts = timeseries(obj.Temp(:,2), obj.Temp(:,1));
    temp = resample(Temp_ts, obj.P(:,1));
    all_temp = [all_temp;temp.data];

    [~, sun] = arrayfun(@(j) obj.sunpos_realworld_v2(j), 1:numel(obj.ti),'UniformOutput', false);
    sun_zths= 100-cell2mat(sun)';
    sun_zths_ts = timeseries(sun_zths, obj.ti);
    sun_zths = resample(sun_zths_ts, obj.P(:,1)); 
    all_sun_zth = [all_sun_zth;sun_zths.data];
    all_pw = [all_pw;obj.P(:,2)];
end
dt = [all_Irr, all_temp, all_sun_zth];
nan_idx = sum(isnan(dt),2);
X=normc(dt(find(nan_idx==0),:));
y=normc(all_pw(find(nan_idx==0)));
% figure;subplot(1,3,1);plot(X(:,1),y);
% subplot(1,3,2);plot(X(:,2),y);
% subplot(1,3,3);plot(X(:,3),y);

% modelFun = @(b,x) b(1).*(x(:,1).^b(2)) + b(3).*(x(:,2).^b(4)) + b(5).*(x(:,3).^b(6)) + b(7);  %a(x^b)
% modelFun = @(b,x) b(1).*x(:,1) + b(2).*x(:,2) + b(3).*x(:,3) + b(4);  %a(x^b)
% beta0 = [0.3,2,0.2,2,0.1,2,0];
% model = fitnlm(X,y,modelFun,beta0);
end

model = fitlm(X,y,'linear','RobustOpts','on');

test_Irr=[];
test_temp=[];
test_sun_zth=[];
test_pw=[];

for d=1:numel(test_days)
    date = char(test_days(d));
    obj = vmlSeq(date);

    Irr_ts = timeseries(obj.Irr(:,2), obj.Irr(:,1));
    Irr = resample(Irr_ts, obj.P(:,1)); 
    test_Irr = [test_Irr;Irr.data];

    Temp_ts = timeseries(obj.Temp(:,2), obj.Temp(:,1));
    temp = resample(Temp_ts, obj.P(:,1));
    test_temp = [test_temp;temp.data];

    [~, sun] = arrayfun(@(j) obj.sunpos_realworld_v2(j), 1:numel(obj.ti),'UniformOutput', false);
    sun_zths= 180-cell2mat(sun)';
    sun_zths_ts = timeseries(sun_zths, obj.ti);
    sun_zths = resample(sun_zths_ts, obj.P(:,1)); 
    test_sun_zth = [test_sun_zth;sun_zths.data];
    test_pw = [test_pw;obj.P(:,2)];
end
dt = [test_Irr, test_temp, test_sun_zth];
nan_idx = sum(isnan(dt),2);
X_t=normc(dt(find(nan_idx==0),:));
y_t=normc(test_pw(find(nan_idx==0)));

est_pw = predict(model,X_t);
est_pw = max(0,est_pw);
figure;
h(1) = plot(est_pw,'r.');
hold on;
h(2) = plot(y_t,'b.');
legend([h(1),h(2)],{'Predict'; 'Observed'});
ylabel('Power');

