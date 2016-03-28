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

% DHI1 = irr1 - cosd(zenith).*DNI;
% DHI2 = irr2 - max(0,cosd(zenith+elev)).*DNI;

cal_irr1 = a1*cosd(zenith).*DNI + b1*DHI;
cal_irr2 = a2*max(0,cosd(angles')).*DNI + b2*DHI;

% ind = 1;
% A=[a1*cosd(zenith(ind)) b1;
%    a2*max(0,cosd(angles(ind))) b2];
% dni_dhi = A\[irr1(ind); irr2(ind)]

figure(107);
subplot(2,1,1);
plot(irr1,cal_irr1,'.');
xlabel('GHI observed');
ylabel('GHI reconstructed');
hold on;
plot([0 1200], [0 1200],'r-');
grid on;
title('Correlation of GHI reconstruction for horizontal sensor');

subplot(2,1,2);
plot(irr2,cal_irr2,'.');
xlabel('GHI observed');
ylabel('GHI reconstructed');
grid on;
hold on;
plot([0 600], [0 600],'r-');
title('Correlation of GHI reconstruction for tilted sensor');
%%%%%%%%%%%%%%%%%%%%%%%%%%5



figure;
subplot(2,1,1);
values = hist3([irr1(:) cal_irr1(:)],[1200 1200]);
newmap = jet;
newmap(1,:) = [1 1 1];
colormap(newmap); 
imagesc(values)
colorbar
% axis equal
axis xy
% plot(irr1,cal_irr1,'.');
xlabel('GHI observed');
ylabel('GHI reconstructed');
hold on;
plot([0 1200], [0 1200],'r-');
grid on;
title('Correlation of GHI reconstruction for horizontal sensor');

subplot(2,1,2);
plot(irr2,cal_irr2,'.');
xlabel('GHI observed');
ylabel('GHI reconstructed');
grid on;
hold on;
plot([0 600], [0 600],'r-');
title('Correlation of GHI reconstruction for tilted sensor');

% values = hist3([data1(:) data2(:)],[51 51]);
% imagesc(values)
% newmap = jet;
% newmap(1,:) = [1 1 1];
% colormap(newmap);    
% colorbar
% axis equal
% axis xy