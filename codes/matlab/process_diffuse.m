img_save_path='';
prj_path='';
proj_path;
   
addpath(prj_path);
conf = [];
conf=local_conf(conf);
load('calc\data_clean.mat', 'data');
model3D=load([conf.datafolder  'Cavriglia_model3D.mat']);
model3D = model3D.model3D;

times = cellfun(@(d) d.time, data);
irr = cellfun(@(d) d.irr(1), data);
irr45 = cellfun(@(d) d.irr(2), data);
diffuse = cellfun(@(d) d.diff_irr, data);
clear_irr = cellfun(@(d) d.clear_irr(2), data);
sun_flag = cellfun(@(d) d.sun_flag, data);
[ClearSkyGHI,ClearSkyDNI,ClearSkyDHI,Zenith,Azimuth] = pvl_clearsky_ineichen(pvl_maketimestruct(times, ...
model3D.UTC),model3D.Location);

az = .7;
elev = 33.5;    
plate_cord = repmat([deg2rad(az),deg2rad(elev),1],[length(times),1]);
[px,py,pz] = sph2cart(plate_cord(:,1),plate_cord(:,2),plate_cord(:,3));
plate_cord = [px,py,pz];
[sx,sy,sz] = sph2cart(deg2rad(Azimuth),deg2rad(90-Zenith),ones(size(times))');
sun_cords = [sx,sy,sz];
nrm=sqrt(sum(abs(cross(sun_cords,plate_cord,2)).^2,2));
angles = atan2d(nrm, dot(sun_cords,plate_cord,2));
effective_DNI = clear_irr'.*max(0,cosd(angles));
other_diffuse = irr45-effective_DNI'.*sun_flag_to_coef(sun_flag);

figure(1); plot(times,diffuse,'b.',times,other_diffuse,'r.--');
ylim = get(gca,'ylim'); ylim(1) = 0;
set(gca,'ylim',ylim);
grid on;
datetickzoom;

for i=1:length(data)
    d=data{i};
    img=imread([img_save_path d.day '__' num2str(d.j) '.jpeg']);
    cloud_shine = double(rgb2gray(img)-150);
    cloud_shine(cloud_shine<0)=0;
    saturation = nansum(nansum(cloud_shine));
    lab_img = colorspace('RGB->LAB',img);        
    saturation2 = sqrt(lab_img(:,:,2).^2 + lab_img(:,:,3).^2);
    saturation2(saturation2<100)=0;
    saturation2 = nansum(nansum(saturation2));
end