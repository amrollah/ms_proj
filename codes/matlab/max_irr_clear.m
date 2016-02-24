close all; 
clear all; 
clc;

%% settings
prj_path='';
proj_path;
   
addpath(prj_path);
conf = [];
conf=local_conf(conf);
model3D=load([conf.datafolder  'Cavriglia_model3D.mat']);
model3D = model3D.model3D;

days = dir(conf.basefolder); % get all content of Cavriglia image folder
days = days(vertcat(days.isdir)); % filter only folders
days = days(8:end-1)';


start_date = '';
values = {};

start_date_found = true;
for d=1:numel(days)
    date = days(d).name;
    if strcmp(start_date, date) || start_date_found
        start_date_found = true;
    else
        continue;
    end
%     try
        disp(date);
        day=[];
        day.date=date;
        folder1 = textscan(date,'%f','delimiter','_');
        dtv = reshape(folder1{1}(1:1:3),1,3);
        times = [];
        times.UTCOffset = [];
        times.year = [];
        times.month = [];
        times.day = [];
        times.hour = [];
        times.minute = [];
        times.second = [];
        for hour=10:14
            for min=0:2:59
                times.UTCOffset(end+1,1) = model3D.UTC;
%                 [Time.year, Time.month, Time.day, Time.hour, Time.minute, Time.second]
                times.year(end+1,1) = dtv(1);
                times.month(end+1,1) = dtv(2);
                times.day(end+1,1) = dtv(3);
                times.hour(end+1,1) = hour;
                times.minute(end+1,1) = min;
                times.second(end+1,1) = 0;
            end
        end
        [ClearSkyGHI,ClearSkyDNI,ClearSkyDHI,Zenith,Azimuth] = pvl_clearsky_ineichen(times,model3D.Location);
        day.max_irr = max(ClearSkyGHI);
%         day.min_irrr = min(ClearSkyGHI);
        values{end+1} = day;
%     catch
%         continue;
%     end
end
save('calc\max_irr.mat', 'values');