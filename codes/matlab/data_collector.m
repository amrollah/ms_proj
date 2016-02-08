close all; 
clear all; 
clc;

%% settings
% prj_path = 'C:\Users\chamsei\Documents\GitHub\ms_proj';
prj_path = 'E:\ms_proj\';
addpath(genpath([prj_path,'\old_codes\vismolib_v5']));
addpath(genpath([prj_path,'\lib']));

cav_img_path = 'E:\ABB\cav\img';
days = dir(cav_img_path); % get all content of Cavriglia image folder
days = days(vertcat(days.isdir)); % filter only folders
days = days(8:end-1)';

step = 1;
counter = 1;

for d=1:numel(days)
    date = days(d).name;
    s = vmlSeq('cavriglia',date);
    jrange = s.find_daylight_range();
    for j=jrange(1):step:jrange(end)
        s.loadframe(j);
        frame.sun_flag = s.calc.seg.clear_sun_flag(j);
        if  frame.sun_flag == 1 || frame.sun_flag == 4
            frame.day=date;
            frame.j=j;
            frame.irr = s.data.Irr(j,2:3);
            frame.diff_irr = s.calc.Irr(j,2);
            frame.clear_irr = s.data.IrrClear(j,2:end);
            frame.zenith = s.data.Zenith(j,2);
            frame.cloud_shine = s.calc.seg.cloud_shine_fact(j);
            frame.clouds = obj.calc.seg.clouds_fact(j);
            frame.power = s.data.P(j,2);
            data{counter}=frame;
            
            counter = counter + 1;
        end
    end
    clear s;
end
save('img_data.mat', data);
