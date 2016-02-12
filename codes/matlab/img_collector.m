close all; 
clear all; 
clc;

%% settings
prj_path = 'C:\Users\chamsei\Documents\GitHub\ms_proj\';
addpath(genpath([prj_path,'\old_codes\vismolib_v5']));
addpath(genpath([prj_path,'\lib']));
addpath(prj_path);

cav_img_path = 'U:\HDRImages_Cavriglia\img\';
img_save_path = 'C:\data\diffuse\';
load('calc\img_data.mat', 'data')
last_day = '';

for i=1:length(data)
    d = data{i};
    if ~strcmp(d.day, last_day)
        last_day=d.day;
        disp(last_day);
        try
            s = vmlSeq('cavriglia',last_day);
        catch
            continue;
        end
    end
    d.time = s.data.ti(d.j);
%     file_name = strcat(s.data.files(d.j),'_Debevec.jpeg');
%     copyfile(char(strcat(cav_img_path,last_day,'\',file_name)),strcat(img_save_path,last_day,'__',num2str(d.j),'.jpeg'));  
    data{i} = d;
end

save('calc\img_data_with_time.mat', 'data');