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
load('calc\img_data.mat', 'data')
last_day = '';
cc = {};
counter = 0;
batch = 11;
for i=28879:length(data)
    try
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
    
%     file_name = strcat(s.data.files(d.j),'_Debevec.jpeg');
%     copyfile(char(strcat(cav_img_path,last_day,'\',file_name)),strcat(img_save_path,last_day,'__',num2str(d.j),'.jpeg'));  
%     if d.diff_irr < .1*d.irr(1) || d.diff_irr > .5*d.clear_irr(2)
%         s.loadframe(d.j);
%         d.sun_flag = s.curseg.clear_sun_flag;
%         d.clouds = s.calc.seg.clouds_fact(d.j);
%         d.diff_irr = d.irr(1)-d.clear_irr(2)*cosd(d.zenith)*s.sunFlagToCoef(d.sun_flag,2);
%     end
        s.loadframe(d.j);        
        frame.day=d.day;    
        frame.j=d.j;
        frame.ind = uint32(find(s.curseg.cc(:)==1));
        cc{end+1} = frame;
        counter = counter + 1;
        if counter==2000
            save(['C:\data\cc_data\cc_data', num2str(batch), '.mat'], 'cc');
            clear cc; cc = {}; counter = 0; batch = batch+1;
        end
    catch
        continue;
    end
end

% save('calc\cc_data.mat', 'cc');