close all; clear all; clc;

%% settings
addpath(genpath('C:\Users\chamsei\Documents\GitHub\ms_proj\old_codes\vismolib'))
vmlconfig_cavriglia;
conf = evalin('base','VMLCONF');
save_figs = true;
show_power = false;

%% Path of data
% folder of McClear export files
mcclear_path = 'C:\data\mcclear\';
% base folder for Cavriglia image and data files
cav_data_path = 'U:\HDRImages_Cavriglia\img\';
local_cav_data_path = 'C:\data\cav\log_data\';

%% read McClear data file (estimated data)
filename_mcclear = strcat(mcclear_path, 'irradiation-clear_sky_cav_2015_till_september_v2.csv');
fid = fopen(filename_mcclear);
mc_data = textscan(fid,'%s%s%f%f%f%f%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s','delimiter',';');
fclose(fid);

%% find all available days for irradiation visualization
days = dir(cav_data_path); % get all content of Cavriglia image folder
days = days(vertcat(days.isdir)); % filter only folders
days = days(8:end-1)'; % skip first 5 folders (first two are . and .., data of next 5 are corrupted) and skip last folder, because log data is not still transfered for current day.
%only consider these days (because they are clear_sky)
days  = {
%     '2015_08_03',
    '2015_08_04',
    '2015_08_26',
    '2015_09_21'
    };

RMSE_mat = {}; % the array to store RMSE values for each day
start_date = '2015_07_26'; %days(1).name; % '2015_08_28'; % if not wish to provide a start_date, just replace the specific date with days(d).name;
end_date = '2015_11_08'; % we don't consider days after this date
start_date_found = false;
day_counter = 1; % counter for days which we decide to compare and store RMSE

%%% loop over all days
for d=1:size(days,1)
% set the date for comparison
if isfield(days, 'name')
    date = days(d).name;
    if strcmp(start_date, date) || start_date_found
    start_date_found = true;
    else
        continue;
    end
else
    date = char(days(d));
end
disp(date);

%% read PV log data file (real data)
log_file = strcat(local_cav_data_path, date, '_data.log');
if ~exist(log_file, 'file')
    filename_log = strcat(cav_data_path, date, '\', date, '_data');
    [status,result] = system(['"C:\Program Files\7-Zip\7z.exe" -y x ' '"' filename_log '.7z"' ' -o' '"' local_cav_data_path '"']); 
end
try
    fid = fopen(log_file);
    log_data = textscan(fid,'%s%s%f%f%f%f%f%f','delimiter',' ');
    fclose(fid);
catch
    % the file does not exist or is corrupted; so we skip this day
    continue;
end

obj = vmlSeq(date,[6 20]);

%% plot complete version of log data in a temp figure just to see if current day is a clear sky or not. 
% We don't continue to next steps for mostly cloudy days.
n = size(log_data{1,1},1);
temp_fig = figure(1000);
%plot(log_data{1,5}, '--g');
plot(obj.Irr(:,2), '--g');
hold on;
plot(obj.ClearSkyRef(:,2), '--b');
hold on;
plot(obj.ClearSkyOrigRef(:,2), '--r');
xtick = 1:500:size(log_data{1,1},1);
set(gca, 'XTick', xtick, 'XTickLabel',log_data{1,2}(xtick));
rotateXLabels(gca, 90)
title(strcat(strrep(date, '_', '/'), ' irradiance'));
xlabel('Time');
ylabel('Irradiance (W/m^2)');

% choice = menu('Do you want to compare this day with McClear data and save the figure as an image.','Yes','No');
% close(1000); % close the temp figure.
% if choice==2 || choice==0
%     continue;
% end;


% visualize clear times and cloudy time
figure;
plot(obj.Irr(obj.clear_times,1),obj.Irr(obj.clear_times,2), '.g');
hold on;
plot(obj.Irr(obj.cloudy_times,1),obj.Irr(obj.cloudy_times,2), '.r');
hold on;
plot(obj.ClearSkyRef(:,1),obj.ClearSkyRef(:,2), '--b');



%% iterate over both files to find matching moments and extract the data
% initalization of params
time_shift = 2;

tim = char(log_data{1,2}(1));
time = tim(1:end-3);
last_time = time;

last_time_shifted = strcat(num2str(str2double(time(1:2))-time_shift), time(3:end));
i = 1;
idx = 1;
counter = 1;
data_est_label = {};
% loop over real data log file (incidents are logged every couple of  seconds)
while i<=n
    avg_real_minute = 0;
    last_i = i;
    while i<=n % aggregate (average) irradiation for every minute
        tim = char(log_data{1,2}(i));
        time = tim(1:end-3);
        avg_real_minute = avg_real_minute + log_data{1,5}(i);
        if ~strcmp(last_time,time) && i~=1
            break;
        end    
        i = i + 1;
    end
    if size(last_time_shifted) < 5 % result of num2str for numbers less then 10 needs to be augmented with a '0' 
        mcclear_tim = strcat('0',last_time_shifted,':00');
    else
        mcclear_tim = strcat(last_time_shifted,':00');
    end
    last_time = time;
    last_time_shifted = strcat(num2str(str2double(time(1:2))-time_shift), time(3:end));
    data_real(counter) = avg_real_minute/(i - last_i);
    % find the matching line in McClear file
    while ~strcmp(mc_data{1,1}(idx), log_data{1,1}(last_i)) || ~strcmp(mc_data{1,2}(idx), mcclear_tim) 
       idx = idx + 1;
    end
    data_est(counter,1) = mc_data{1,5}(idx)*60; % direct beam irradiation
    data_est(counter,2) = mc_data{1,6}(idx)*60; % diffuse irradiation
    data_est_label{counter} = last_time;
    counter = counter + 1;
end
adj_data_est = adjust_mcclear(data_est(:,1) + data_est(:,2));
RMSE_mat{day_counter,1} = date;

% generate reference clear sky irradiation
LinkeTurbidity = prep_LinkeTurbidity(obj);
[ClearSkyGHI, ClearSkyDNI, ClearSkyDHI] = arrayfun(@(tm) pvl_clearsky_ineichen(time_struct(obj,tm), obj.calib.model3D.Location,'LinkeTurbidityInput', LinkeTurbidity), data_est_label(1,:));

choice = menu('Continue?','Yes');
continue;
%% plot real log data and mcclear estimated data with reference from PVlib for current day
f = figure(d);
h(1) = plot((1:counter-1),data_real(1:counter-1), '--b');
hold on;
ylim([0, max(max(data_real), max(adj_data_est)) + 30]);
h(2) = plot((1:counter-1),adj_data_est(1:counter-1), '--r');
h(3) = plot((1:counter-1),adjust_reference(ClearSkyGHI), '--g');
%h(3) = plot((1:counter-1),adjust_mcclear(data_est(1:counter-1,2)), '--g');
legend([h(1),h(2),h(3)],{'Real data'; 'McClear'; 'Reference'}); 
xtick2 = 1:100:size(data_est_label,2);
set(gca, 'XTick', xtick2, 'XTickLabel',data_est_label(1,xtick2));
rotateXLabels(gca, 90)
title(strcat(strrep(date, '_', '/'), ' irradiance'));
xlabel('Time');
ylabel('Irradiance (W/m^2)');

hold off;

if show_power
%% ploting the power output
    figure(2000);
    plot(obj.P(:,2), '--b');
    xtick = 1:300:size(obj.P(:,2),1);
    set(gca, 'XTick', xtick, 'XTickLabel', arrayfun(@(tm) time_convertor(datevec(tm)), obj.P(xtick,1), 'UniformOutput', false));
    rotateXLabels(gca, 90)
    title(strcat(strrep(date, '_', '/'), ' power output'));
    xlabel('Time');
    ylabel('Power (W)');
    choice = menu('Continue?','Yes');
    close(2000); % close the temp figure.
end

begin_t = 1; %100; %beginging of time index
end_t = numel(data_real); %300; %end of time index
% calc the RMSE error for estimated data
rmse=sqrt(mean((data_real(begin_t:end_t)'-adj_data_est(begin_t:end_t,1)).^2));
rmse_percent = rmse/mean(data_real(begin_t:end_t))*100;
RMSE_mat{day_counter,2} = rmse_percent;
day_counter = day_counter + 1;
fprintf('RMSE: %f\n', rmse_percent);

if (save_figs)
    saveas(f, strcat(local_cav_data_path, date, '.jpg'));
end
if strcmp(end_date, date)
     break;
end
end
save('RMSE_mat.mat', 'RMSE_mat');