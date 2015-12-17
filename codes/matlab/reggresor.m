close all; 
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

test_days={'2015_08_04',};
load_data = true;

if load_data
%     load('lnModel.mat');
    load('Data.mat');
else
all_Irr=[];
all_temp=[];
all_sun_zth=[];
all_pw=[];
for d=1:numel(days)
    date = char(days(d));
    obj = vmlSeq(date);

    Irr = max(0,interp1(obj.Irr(:,1),obj.Irr(:,2),obj.P(:,1),'linear'));
    all_Irr = [all_Irr;Irr];
    temp = interp1(obj.Temp(:,1),obj.Temp(:,2),obj.P(:,1),'linear');
    all_temp = [all_temp;temp];

    [~, sun] = arrayfun(@(j) obj.sunpos_realworld_v3(j), obj.Irr(:,1),'UniformOutput', false);
    sun_zths= 180-cell2mat(sun)';
    sun_zths = interp1(obj.Irr(:,1),sun_zths,obj.P(:,1),'linear');
    all_sun_zth = [all_sun_zth;sun_zths];
    all_pw = [all_pw;obj.P(:,2)];
    figure;
    plot(Irr,'k.');
    hold on;
    plot(temp,'b.');
    hold on;
    plot(sun_zths,'g.');
    pause(1);
end
dt = [all_Irr, all_temp, all_sun_zth];
nan_idx = sum(isnan(dt),2);
X=dt(find(nan_idx==0),:);
y=all_pw(find(nan_idx==0));
% figure;subplot(1,3,1);plot(X(:,1),y);
% subplot(1,3,2);plot(X(:,2),y);
% subplot(1,3,3);plot(X(:,3),y);

% modelFun = @(b,x) b(1).*(x(:,1).^b(2)) + b(3).*(x(:,2).^b(4)) + b(5).*(x(:,3).^b(6)) + b(7);  %a(x^b)
% modelFun = @(b,x) b(1).*x(:,1) + b(2).*x(:,2) + b(3).*x(:,3) + b(4);  %a(x^b)
% beta0 = [0.3,2,0.2,2,0.1,2,0];
% model = fitnlm(X,y,modelFun,beta0);


test_Irr=[];
test_temp=[];
test_sun_zth=[];
test_pw=[];

for d=1:numel(test_days)
    date = char(test_days(d));
    obj = vmlSeq(date);

    Irr = max(0,interp1(obj.Irr(:,1),obj.Irr(:,2),obj.P(:,1),'linear'));
    test_Irr = [test_Irr;Irr];

    temp = interp1(obj.Temp(:,1),obj.Temp(:,2),obj.P(:,1),'linear');
    test_temp = [test_temp;temp];

    [~, sun] = arrayfun(@(j) obj.sunpos_realworld_v3(j), obj.Irr(:,1),'UniformOutput', false);
    sun_zths= 180-cell2mat(sun)';
    sun_zths = interp1(obj.Irr(:,1),sun_zths,obj.P(:,1),'spline');
    test_sun_zth = [test_sun_zth;sun_zths];
    test_pw = [test_pw;obj.P(:,2)];
end
dt = [test_Irr, test_temp, test_sun_zth];
nan_idx = sum(isnan(dt),2);
X_t=dt(find(nan_idx==0),:);
y_t=test_pw(find(nan_idx==0));

save('Data.mat', 'X','y','X_t','y_t');
end

X=normc(X);
Y=normc(Y);
X = [X(:,[1,3]), sqrt(X(:,2))];
X_t = [X_t(:,[1,3]), sqrt(X_t(:,2))];
X_t=normc(normc);
y_t=normc(y_t);
model = fitlm(X,y,'linear','RobustOpts','on');

est_pw = predict(model,X_t);
est_pw = max(0,est_pw);
figure;
h(1) = plot(est_pw,'r.');
hold on;
h(2) = plot(y_t,'b.');
legend([h(1),h(2)],{'Predict'; 'Observed'});
ylabel('Power');

est_pw = predict(model,X);
est_pw = max(0,est_pw);
figure;
h(1) = plot(est_pw,'r.');
hold on;
h(2) = plot(y,'b.');
legend([h(1),h(2)],{'Predict'; 'Observed'});
ylabel('Power');
hold on;
plot(X(:,1),'k.');
hold on;
plot(X(:,2),'c.');
hold on;
% plot(X(:,3),'y.');
% hold on;
% plot(X(:,4),'g.');
