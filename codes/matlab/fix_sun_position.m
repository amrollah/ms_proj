close all;
clear all;
clc;
img_save_path='';
prj_path='';
proj_path;

addpath(prj_path);
s = vml('cavriglia','2015_12_11',[],true);
% R=s.kabsch([],[500,800,1000,1300,1600,2000,2300,2500,2800,3000,3300,3500,3700]);
load('calc\clean_data_with_8cc_nan_corrected2.mat', 'data');
p_w = 4;
R = 80;
    %29282
for i=1:length(data)
        d=data{29282};
%         if d.sun_flag==1 || ~strcmp('2015_12_', d.day(1:8))
%             continue;
%         end
        disp(i);
        tm = datevec(d.time);
        img_file = [img_save_path d.day '__' num2str(d.j) '.jpeg'];
        I = s.get_image(img_file);
        sun_pos = s.sunpos_im(d.time);
        cx=floor(sun_pos(1)); cy=floor(sun_pos(2));
        xgrid = max(cx-R,1):p_w:min(cx+R,size(I,1));
        ygrid = max(cy-R,1):p_w:min(cy+R,size(I,2));
        [X,Y] = meshgrid(xgrid,ygrid);
        area = imcrop(I,[ygrid(1),xgrid(1),xgrid(end)-xgrid(1),xgrid(end)-xgrid(1)]);
        h=figure(1); imshow(I); hold on; plot(cy,cx,'ro','MarkerSize',10);
        maxfig(h,1);
        pause(2);
end
