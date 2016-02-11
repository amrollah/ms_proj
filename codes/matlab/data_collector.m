close all; 
% clear all; 
clc;

%% settings
% prj_path = 'C:\Users\chamsei\Documents\GitHub\ms_proj';
prj_path = 'C:\Users\chamsei\Documents\GitHub\ms_proj';
addpath(genpath([prj_path,'\old_codes\vismolib_v5']));
addpath(genpath([prj_path,'\lib']));

cav_img_path = 'U:\HDRImages_Cavriglia\img';
days = dir(cav_img_path); % get all content of Cavriglia image folder
days = days(vertcat(days.isdir)); % filter only folders
days = days(8:end-1)';

start_date = '2015_10_30';
init_step = 2;
std_thres = 0.15;

start_date_found = false;
for d=1:numel(days)
    date = days(d).name;
    if strcmp(start_date, date) || start_date_found
        start_date_found = true;
    else
        continue;
    end
    try
        s = vmlSeq('cavriglia',date);
    catch
        continue;
    end
    jrange = s.find_daylight_range();
    step = init_step;
    j=jrange(1);
    while j<jrange(end) && j<size(s.data.Irr,1)
        tid=s.getClearId(j);
        if abs(s.data.Irr(j,2) - s.data.IrrClear(tid,2)) > 0.12*s.data.IrrClear(tid,2) && s.data.IrrClear(tid,2) > 100
            s.loadframe(j);
            frame.sun_flag = s.calc.seg.clear_sun_flag(j);
            if  frame.sun_flag == 1 || frame.sun_flag == 4

                frame.day=date;
                frame.j=j;
                frame.irr = s.data.Irr(j,2:3);
                frame.diff_irr = s.calc.Irr(j,2);
                frame.clear_irr = s.data.IrrClear(tid,2:end);
                frame.zenith = s.data.Zenith(tid,2);
                frame.cloud_shine = s.calc.seg.cloud_shine_fact(j);
                frame.clouds = s.calc.seg.clouds_fact(j);
                frame.power = s.getP(s.data.ti(j));
                data{end+1}=frame;
            end
        end
        dt = std(s.data.Irr(max(1,j-60):j,2))/s.data.Irr(j,2);
        if dt < std_thres      
            step = min(step*2,30);
        else
            step = init_step;
        end
        j = j+ step;
    end
    clear s;
end
save('img_data.mat', 'data');
