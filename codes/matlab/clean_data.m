clear all;
close all;
clc;
img_save_path='';
prj_path='';
proj_path; 

% az = .7;
% elev = 33.5;
%       
conf = [];
conf=local_conf(conf);
load('calc\clean_data_with_8cc_nan_corrected2.mat', 'data');
% model3D=load([conf.datafolder  'Cavriglia_model3D.mat']);
% model3D = model3D.model3D;

% y = cellfun(@(d) d.corr_tilt_diff, data);

% times = cellfun(@(d) d.time, data);
%     irr = cellfun(@(d) d.irr(1), data);
% irr45 = cellfun(@(d) d.irr(2), data);
% diffuse = cellfun(@(d) d.diff_irr, data);
% tilted_diffuse = cellfun(@(d) d.tilt_diff, data);
% clear_irr = cellfun(@(d) d.clear_irr(2), data);
% sun_flag = cellfun(@(d) d.sun_flag, data);
% [ClearSkyGHI,ClearSkyDNI,ClearSkyDHI,Zenith,Azimuth] = pvl_clearsky_ineichen(pvl_maketimestruct(times, ...
% model3D.UTC),model3D.Location);
% 
% plate_cord = repmat([deg2rad(az),deg2rad(elev),1],[length(times),1]);
% [px,py,pz] = sph2cart(plate_cord(:,1),plate_cord(:,2),plate_cord(:,3));
% plate_cord = [px,py,pz];
% [sx,sy,sz] = sph2cart(deg2rad(Azimuth),deg2rad(90-Zenith),ones(size(times))');
% sun_cords = [sx,sy,sz];
% nrm=sqrt(sum(abs(cross(sun_cords,plate_cord,2)).^2,2));
% angles = atan2d(nrm, dot(sun_cords,plate_cord,2));
% effective_DNI = ClearSkyDNI.*max(0,cosd(angles));
% tilted_diffuse = irr45-effective_DNI'.*sun_flag_to_coef(sun_flag);
% 
% rel_err = 100*(diffuse-tilted_diffuse)./diffuse;

% diffuse2 = medfilt1(diffuse);
% correction on tilted diffuse
% corrected_tilted_diffuse = tilted_diffuse*1.3;
% load('calc\max_irr.mat', 'values');
cl_data = {};
day='';
counter = 0;
for i=1:length(data)
    d = data{i};
%     img = imread([img_save_path d.day '__' num2str(d.j) '.jpeg']);
%     red_m=mean(mean(img(:,:,1)));
%     if red_m>100
%        disp([num2str(i), '     ', num2str(red_m)]);
%     end
%     if red_m<150 && d.diff_irr > .05*d.irr(1) && d.diff_irr < .6*d.clear_irr(2)
%         figure(1);
%         imshow(img);
%     if ~(strcmp(d.day,'2015_09_15') || strcmp(d.day,'2015_09_16') || strcmp(d.day,'2015_09_17') || strcmp(d.day,'2015_09_18')...
%             || strcmp(d.day,'2015_10_27') || strcmp(d.day,'2015_10_28') || strcmp(d.day,'2015_10_29') || strcmp(d.day,'2015_10_30')...
%             || strcmp(d.day,'2015_10_31') || strcmp(d.day,'2015_11_01') || strcmp(d.day,'2015_11_02') || strcmp(d.day,'2015_11_03')...
%             || strcmp(d.day,'2015_11_04') || strcmp(d.day,'2015_11_05') || strcmp(d.day,'2015_11_06') || strcmp(d.day,'2015_11_07')...
%             || strcmp(d.day,'2015_11_08') || strcmp(d.day,'2015_11_09') || strcmp(d.day,'2016_01_02'))
%         cl_data{end+1}=d;
%     end
%     if ~strcmp(d.day,day)
%         day = d.day;
%         disp(day);
%         for j=1:length(values)
%             if strcmp(values{j}.date, day)
%                 max_irr= values{j}.max_irr;
%             end
%         end
%     end
%     d.tilt_diff = tilted_diffuse(i);
%     d.azimuth = Azimuth(i);
%     if d.clear_irr(1)/max_irr > .25
%     if abs(rel_err(i)) < 160
%        t = datevec(d.time);
%        disp(t([4,5]));
%        counter = counter + 1;
%        cl_data{end+1}=d;
%     end
% d.corr_tilt_diff = corrected_tilted_diffuse(i);
% d.diff_median = diffuse2(i);
% if d.corr_tilt_diff>10
%     cl_data{end+1}=d;
% end
if ~isnan(mean(d.cc_fact))
    cl_data{end+1}=d;
end
end
data=cl_data;
%44519
save('calc\clean_data_with_8cc_nan_corrected3.mat', 'data');