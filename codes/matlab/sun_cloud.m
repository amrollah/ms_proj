close all; 
clear all; 
% clc;

%% settings
addpath(genpath('C:\Users\chamsei\Documents\GitHub\ms_proj\old_codes\vismolib'))
addpath(genpath('C:\Users\chamsei\Documents\GitHub\ms_proj\old_codes\solvers'));
vmlconfig_cavriglia;
conf = evalin('base','VMLCONF');


%% Path of data
% base folder for Cavriglia image and data files
cav_data_path = 'U:\HDRImages_Cavriglia\img\';
local_cav_data_path = 'C:\data\cav\suns\';

suns = dir('C:\data\cav\sun\');
suns = suns(~vertcat(suns.isdir));
start_date = '2015_08_03'; 
start_date_found = false;

dstr='yyyy_mm_dd_HH_MM_SS';
last_date = '';
% loop
for d=1:numel(suns)
    if isfield(suns, 'name')
        orig_im = suns(d).name;
        im = strrep(orig_im,'_Debevec.jpeg','');
        date = im(1:10);
        if strcmp(start_date, date) || start_date_found
            start_date_found = true;
        else
            continue;
        end
    else
        im = char(suns(d));
    end
    disp(im);
    
    if ~strcmp(last_date,date)
        s = vmlSeq(date);
        last_date = date;
    end
    tm = datenum(im,dstr);
    [~,j] = min(abs(s.ti-tm));
    s.showframe(j);
    s.loadframe(j);
    sat_sun_pos = obj.detect_saturated_sun(j);
    sun_pos = s.sun_pos([2 3],j);
    f=s.xcur.r2b>=s.xcur.thres.vmfi;
    g=imresize(f,s.sz);
    mask = s.skymask_wo_sun(j);
    mask(isnan(mask)) = false;
    g(~mask) = 0;
    P = s.plant_projection_on_image(j,1000);
    [x,y]=find(g>-1);
%     figure(110);imshow(g)
    g2 = inpolygon(x,y,P(1,:),P(2,:));
    g2 = reshape(g2,size(g));
    shadow = g2&g;
    figure(130);imshow(shadow)

    
    pause(1);
end


