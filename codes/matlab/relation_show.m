function relation_show()
close all; 
% clear all; 
% clc;
img_save_path='';
prj_path='';
proj_path;
addpath(prj_path);
if ~exist('data','var')
    load('calc\clean_data_with_8cc_nan_corrected3.mat', 'data');
    DHI = cellfun(@(d) d.corr_tilt_diff, data);
    sun_flag = cellfun(@(d) d.sun_flag, data);
    sun_flag=normalize(sun_flag,0,1);
    clear_DNI = cellfun(@(d) d.clear_irr(2), data);
    DNI = clear_DNI .* sun_flag;
    zenith = cellfun(@(d) d.zenith, data);
    azimuth = cellfun(@(d) d.azimuth, data);
    irr1 = cellfun(@(d) d.irr(1), data);
    irr2 = cellfun(@(d) d.irr(2), data);
end
a1=.9;
b1=1.14;
a2=1;
b2=.78;
elev=33.5;

az = 12;   
plate_cord = repmat([deg2rad(az),deg2rad(elev),1],[length(irr1),1]);
[px,py,pz] = sph2cart(plate_cord(:,1),plate_cord(:,2),plate_cord(:,3));
plate_cord = [px,py,pz];
[sx,sy,sz] = sph2cart(deg2rad(azimuth),deg2rad(90-zenith),ones(size(irr1)));
sun_cords = [sx',sy',sz'];
nrm=sqrt(sum(abs(cross(sun_cords,plate_cord,2)).^2,2));
angles = atan2d(nrm, dot(sun_cords,plate_cord,2));
% effective_DNI = clear_irr'.*max(0,cosd(angles));

DHI1 = irr1 - cosd(zenith).*DNI;
DHI2 = irr2 - max(0,cosd(zenith+elev)).*DNI;

cal_irr1 = a1*cosd(zenith).*DNI + b1*DHI;
cal_irr2 = a2*max(0,cosd(angles')).*DNI + b2*DHI;
figure(107);
subplot(2,2,2);
plot(irr1,cal_irr1,'.');
xlabel('irr horiz');
ylabel('calced_irr_horiz');
hold on;
plot([0 1200], [0 1200],'r-');
grid on;
title('irr1 vs calc_irr1');
subplot(2,2,4);
plot(irr2,cal_irr2,'.');
xlabel('irr tilted');
ylabel('calced_irr tilted');
grid on;
hold on;
plot([0 600], [0 600],'r-');
title('irr2 vs calc_irr2');

show_upd(1);
datacursormode on;
dcmobj = datacursormode(gcf); 
set(dcmobj,'UpdateFcn',@show_upd,'DisplayStyle','window');
function out = show_upd(i, event_obj)
    if nargin==2, i=event_obj.DataIndex; end
    d = data{i};
    figure(107); 
    subplot(2,2,[1,3]); 
    imshow([img_save_path d.day '__' num2str(d.j) '.jpeg']);
    drawnow;
    if nargin<2, out=[]; 
    else
        out = {datestr(d.time,'yyyy-mm-dd HH:MM:SS'), ['i: ', num2str(i)], ['clouds: ',num2str(d.clouds),'%'], ['sat_fact: ', num2str(d.sat_fact)],['sun_flag: ',num2str(d.sun_flag)],['irr1: ', num2str(d.irr(1))],['irr2: ', num2str(d.irr(2))],['diffuse: ', num2str(d.diff_irr)], ['tilted_diffuse: ', num2str(d.corr_tilt_diff)], ['clear_diffuse: ', num2str(d.clear_irr(3))],['DNI: ', num2str(d.clear_irr(2))],['zenith: ',num2str(d.zenith)]};
    end
end
end