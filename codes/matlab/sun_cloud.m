% close all; 
clear all; 
% clc;

%% settings
addpath(genpath('C:\Users\chamsei\Documents\GitHub\ms_proj\old_codes\vismolib'))
addpath(genpath('C:\Users\chamsei\Documents\GitHub\ms_proj\old_codes\solvers'));
vmlconfig_cavriglia;
conf = evalin('base','VMLCONF');
show_plant=false;

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
sun_patches = {};
% loop
f=figure(1);
%  maxfig(f,1);
% set(f,'units','normalized','outerposition',[0 0 1 1]);

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
%     s.showframe(j);
    s.loadframe(j);
    sat_sun_pos = s.detect_saturated_sun(j);
    sun_pos = s.sun_pos([2 3],j);
    if norm(sat_sun_pos-sun_pos)<30
        sun_pos = sat_sun_pos;
    end
   
    cx=sun_pos(1); cy=sun_pos(2); R=20;   
%     cc=1:s.sz(2); 
%     rr=(1:s.sz(1))';
%     f=@(xx,yy) (xx-cx1).^2+(yy-cy1).^2 <=R^2; 
%     C=bsxfun(f,rr,cc);
%     plot(256*C(1),256*C(2),'.b');
%     Inew = im.*repmat(unit8(C),[1,1,3]);

    im = s.imread(j);
    sun_patch = im(cx-R:cx+R,cy-R:cy+R,:); 
%     figure;
%     subplot(14,15,d);
    imshow(sun_patch);
    set(f,'position',get(0,'screensize'));
    sun_rect = getrect(f);
    sun_patches{1,d} = sun_patch;
    sun_patches{2,d} = sun_rect;
    
    if show_plant
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
        figure(130);imshow(shadow);
    end
    pause(1);
end


