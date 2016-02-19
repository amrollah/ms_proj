close all; 
clear all; 
clc;

%% settings
prj_path = 'C:\Users\chamsei\Documents\GitHub\ms_proj\';
addpath(genpath([prj_path,'\old_codes\vismolib_v6']));
addpath(genpath([prj_path,'\lib']));
addpath(prj_path);
cav_img_path = 'U:\HDRImages_Cavriglia\img\';

img_save_path='';
prj_path='';
proj_path;
load('calc\img_data4.mat', 'data')
last_day = '';

for i=57996:length(data)%55855
    d = data{i};
    if ~strcmp(d.day, last_day)
        last_day=d.day;
        disp(last_day);
        try
            s = vml('cavriglia',last_day);
        catch
            continue;
        end
    end
    file_name = strcat(s.data.files(d.j),'_Debevec.jpeg');
    copyfile(char(strcat(cav_img_path,last_day,'\',file_name)),strcat(img_save_path,last_day,'__',num2str(d.j),'.jpeg'));  
%     if d.diff_irr < .1*d.irr(1) || d.diff_irr > .5*d.clear_irr(2)
%         s.loadframe(d.j);
%         d.sun_flag = s.curseg.clear_sun_flag;
%         d.clouds = s.calc.seg.clouds_fact(d.j);
%         d.diff_irr = d.irr(1)-d.clear_irr(2)*cosd(d.zenith)*s.sunFlagToCoef(d.sun_flag,2);
%     end
%     data{i} = d;
end

% save('calc\img_data2.mat', 'data');