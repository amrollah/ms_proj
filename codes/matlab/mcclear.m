close all; clear all; clc;

%% Path of data
% folder of McClear export files
mcclear_path = 'E:\ABB\';
% base folder for Cavriglia image and data files
cav_data_path = 'U:\HDRImages_Cavriglia\img\';
local_cav_data_path = 'E:\ABB\cav\img\';

%% set the date for comparison
date = '2015_08_03';

%% read McClear data file (estimated data)
filename_mcclear = strcat(mcclear_path, 'McClear_July_to_Feb.csv');
fid = fopen(filename_mcclear);
mc_data = textscan(fid,'%s%s%s%s%f%f%f%f%f','delimiter',';');
fclose(fid);

%% read PV log data file (real data)
log_file = strcat(local_cav_data_path, date, '\', date, '_data.log');
% if ~exist(log_file, 'file')
%     filename_log = strcat(cav_data_path, date, '\', date, '_data');
%     [status,result] = system(['"C:\Program Files\7-Zip\7z.exe" -y x ' '"' filename_log '.7z"' ' -o' '"' local_cav_data_path '"']); 
% end
fid = fopen(log_file);
log_data = textscan(fid,'%s%s%f%f%f%f%f%f','delimiter',' ');
fclose(fid);

%% iterate over both files to find matching moments and extract the data
% initalization of params
time_shift = 2;
n = size(log_data{1,1},1);
tim = char(log_data{1,2}(1));
time = tim(1:end-2);
last_time = time;
last_time_shifted = strcat(num2str(str2double(time(1:2))-time_shift), time(3:end));
i = 1;
idx = 1;
counter = 1;
% loop over real data log file (incidents are logged every couple of
% seconds)
while i<=n
    avg_real_minute = 0;
    last_i = i;
    while i<=n % aggregate (average) irradiation for every minute
        tim = char(log_data{1,2}(i));
        time = tim(1:end-2);
        avg_real_minute = avg_real_minute + log_data{1,5}(i);
        if ~strcmp(last_time,time) && i~=1
            break;
        end    
        i = i + 1;
    end
    if size(last_time_shifted) < 6
        mcclear_tim = strcat('0',last_time_shifted,'00');
    else
        mcclear_tim = strcat(last_time_shifted,'00');
    end
    disp(last_time_shifted);
    last_time = time;
    last_time_shifted = strcat(num2str(str2double(time(1:2))-time_shift), time(3:end));
    data_real(counter) = avg_real_minute/(i - last_i);
    % find the matching line in McClear file
    while ~strcmp(mc_data{1,1}(idx), log_data{1,1}(last_i)) || ~strcmp(mc_data{1,2}(idx), mcclear_tim) 
       idx = idx + 1;
    end
    data_est(counter,1) = mc_data{1,7}(idx)*60; % direct beam irradiation
    data_est(counter,2) = mc_data{1,8}(idx)*60; % diffuse irradiation
    data_est(counter,3) = mc_data{1,6}(idx)*60; % total irradiation
    counter = counter + 1;
end


%% plot both real and estimated data for comparison
figure;
plot((1:counter-1),adjust_mcclear(data_est(:,1)+data_est(:,2)), '--r');
hold on;
plot((1:counter-1),data_real, '--b');

begin_t = 100; %beginging of time index
end_t = 300; %end of time index
% calc the RMSE error for estimated data
rmse=sqrt(mean((data_real(begin_t:end_t)'-adjust_mcclear(data_est(begin_t:end_t,1)+data_est(begin_t:end_t,2))).^2));
rmse_percent = rmse/mean(data_real(begin_t:end_t))*100;
 