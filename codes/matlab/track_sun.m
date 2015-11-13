close all;clear all;clc;
addpath(genpath('C:\Users\chamsei\Documents\GitHub\ms_proj\old_codes\vismolib'))
vmlconfig_cavriglia;
conf = evalin('base','VMLCONF');
%load([conf.datafolder 'calib.mat']);
load([conf.datafolder conf.calibration{2}]);

obj = vmlSeq('2015_08_04',[8 18]);
visualize = false;

total_rmse = 0;
counter = 0;
all_pos = [];
for j=1:80:length(obj.ti)
    counter=counter+1;
    [pos, im_pos, pos_m] = sun_position_v2(obj,j,obj.ext_calib.R);
    %im_pos = obj.detect_saturated_sun(j);
    rmse = sqrt(sum((im_pos-pos).^2));
      
    if (visualize)
        figure(counter);
        obj.showframe(j);
        hold on;
        plot(pos(2),pos(1),'o','markersize',10);
        plot(im_pos(2),im_pos(1),'ro','markersize',10);
    end
    all_pos(counter,:) = [pos',im_pos'];
    sign_dist=sign(pos_m-pos);
    fprintf('%d: RMSE:%.2f  dist:%f\n', counter, rmse, sign_dist(1)*pdist([pos,pos_m]','euclidean'));
    total_rmse = total_rmse + rmse;
%     im = cloud_detector( obj, j, pos);
%     imshow(im);
     %pause(1);
end
%save('all_sun_pos.mat', 'all_pos');
fprintf('Avg. RMSE:%.2f  \n', total_rmse/counter);



    