close all;clear all;clc;
addpath(genpath('C:\Users\chamsei\Documents\GitHub\ms_proj\old_codes\vismolib'))
vmlconfig_cavriglia;
conf = evalin('base','VMLCONF');
load([conf.datafolder 'calib.mat']);
obj = vmlSeq('2015_09_27',[8 18]);
visualize = true;

pos_m = sunpos_midday(obj);
x = obj.imread(100);
mid = size(x)/2;
total_rmse = 0;
counter = 0;
thresold = 15;
for j=1:900:length(obj.ti)
    counter=counter+1;
    Q = obj.sunpos_realworld(j);
    pos = obj.camworld2im(R'*Q);
    pos = sun_pos_adjuster( pos, pos_m );
    sign_dist=sign(pos_m-pos);
    im_pos = obj.detect_saturated_sun(j);
    rmse = sqrt(sum((im_pos-pos).^2));
    if rmse<thresold
        pos = im_pos;
        rmse = 0;
    end
    cloud_detector(obj,j);
    obj.opt_cur_thresv();
    if (visualize)
        figure(counter);
        obj.showframe(j);
        hold on;
        plot(pos(2),pos(1),'o','markersize',10);
        plot(im_pos(2),im_pos(1),'ro','markersize',10);
    end
    
    fprintf('%d: RMSE:%.2f  dist:%f\n', counter, rmse, sign_dist(1)*pdist([pos,pos_m]','euclidean'));
    total_rmse = total_rmse + rmse;
    %pause(0.1);
end

fprintf('Avg. RMSE:%.2f  \n', total_rmse/counter);



    