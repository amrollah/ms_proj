close all;
clear all;
clc;
img_save_path='';
prj_path='';
proj_path;

addpath(prj_path);
s = vml('cavriglia','2015_07_30',[],true);
load('calc\data_clean.mat', 'data');
p_w = 4;
R = 60;
    
for i=1:length(data)
        disp(i);
        d=data{i};
        tm = datevec(d.time);
%         if d.sun_flag==4 || tm(4)<12 || tm(4)>15
%             continue;
%         end
        img_file = [img_save_path d.day '__' num2str(d.j) '.jpeg'];
        I = s.get_image(img_file);
        sun_pos = s.sunpos_im(d.time);
        cx=floor(sun_pos(1)); cy=floor(sun_pos(2));
        xgrid = max(cx-R,1):p_w:min(cx+R,size(I,1));
        ygrid = max(cy-R,1):p_w:min(cy+R,size(I,2));
        [X,Y] = meshgrid(xgrid,ygrid);
        area = imcrop(I,[ygrid(1),xgrid(1),xgrid(end)-xgrid(1),xgrid(end)-xgrid(1)]);
        g_area = rgb2gray(area);
        
        sat_contour = g_area>220;
        sat_contour = bwareaopen(sat_contour,10);
        dw=zeros(size(g_area));
        u=zeros(size(g_area));
        r=zeros(size(g_area));
        l=zeros(size(g_area));
        inside_sat = false(size(g_area));
        for xi=1:size(g_area,1)
            dw(xi,:) = sum(sat_contour(xi:end,:),1);
            u(xi,:) = sum(sat_contour(1:xi,:),1);
        end
        for yi=1:size(g_area,2)
            r(:,yi) = sum(sat_contour(:,yi:end),2);
            l(:,yi) = sum(sat_contour(:,1:yi),2);
        end
        for xi=1:size(g_area,1)
            for yi=1:size(g_area,2)
                if ~sat_contour(xi,yi) && sum(r(max(1,xi-1):min(xi+1,size(g_area,1)),yi))>0 && sum(l(max(1,xi-1):min(xi+1,size(g_area,1)),yi))>0 && sum(dw(xi,max(yi-1,1):min(yi+1,size(g_area,2))))>0 && sum(u(xi,max(yi-1,1):min(yi+1,size(g_area,2))))>0
                    inside_sat(xi,yi) = true;
                end
            end
        end
        inside_sat = bwareaopen(inside_sat,10);
        d.sat_fact = min(.25,sum(sum(inside_sat))/numel(g_area));
        ndata{i} = d;
end
data=ndata;
save('calc\data_with_sat_fact.mat', 'data');