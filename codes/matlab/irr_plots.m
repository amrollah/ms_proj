close all; 
clear all; 
clc;

%% settings
prj_path = 'E:\ms_proj';
% prj_path = 'E:\ms_proj\';
addpath(genpath([prj_path,'\old_codes\vismolib_v5']));
addpath(genpath([prj_path,'\lib']));
addpath(prj_path);

cav_img_path = 'E:\ABB\cav\img\';
days = dir(cav_img_path); % get all content of Cavriglia image folder
days = days(vertcat(days.isdir)); % filter only folders
days = days(8:end-1)';


start_date = '2015_08_03';

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
        s = vml('cavriglia',date);
        if size(s.data.Irr,1)>0
            figure(1); s.plot45irr;
            title(strcat(strrep(date, '_', '/'), ' irradiance'));
            saveas(gcf,[s.conf.tmpfolder, 'irr_plots\', date,'.png']);
        end
%     catch
%         continue;
%     end
end