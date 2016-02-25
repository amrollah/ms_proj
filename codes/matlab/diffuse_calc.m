function diffuse_calc(az,elev)
    if nargin<2
      az = .7;
      elev = 33.5;
    end
    img_save_path='';
    prj_path='';
    proj_path;
    
    addpath(prj_path);
    conf = [];
    conf=local_conf(conf);
    load('calc\data_clean.mat', 'data');
%     model3D=load([conf.datafolder  'Cavriglia_model3D.mat']);
%     model3D = model3D.model3D;

    times = cellfun(@(d) d.time, data);
%     irr = cellfun(@(d) d.irr(1), data);
%     irr45 = cellfun(@(d) d.irr(2), data);
    diffuse = cellfun(@(d) d.diff_irr, data);
    diffuse_median = cellfun(@(d) d.diff_median, data);
    tilted_diffuse = cellfun(@(d) d.tilt_diff, data);
    corrected_tilted_diffuse = cellfun(@(d) d.corr_tilt_diff, data);
%     clear_irr = cellfun(@(d) d.clear_irr(2), data);
%     sun_flag = cellfun(@(d) d.sun_flag, data);
%     [ClearSkyGHI,ClearSkyDNI,ClearSkyDHI,Zenith,Azimuth] = pvl_clearsky_ineichen(pvl_maketimestruct(times, ...
%     model3D.UTC),model3D.Location);
%         
%     plate_cord = repmat([deg2rad(az),deg2rad(elev),1],[length(times),1]);
%     [px,py,pz] = sph2cart(plate_cord(:,1),plate_cord(:,2),plate_cord(:,3));
%     plate_cord = [px,py,pz];
%     [sx,sy,sz] = sph2cart(deg2rad(Azimuth),deg2rad(90-Zenith),ones(size(times))');
%     sun_cords = [sx,sy,sz];
%     nrm=sqrt(sum(abs(cross(sun_cords,plate_cord,2)).^2,2));
%     angles = atan2d(nrm, dot(sun_cords,plate_cord,2));
%     effective_DNI = ClearSkyDNI.*max(0,cosd(angles));
%     tilted_diffuse = irr45-effective_DNI'.*sun_flag_to_coef(sun_flag);
%     
    figure; plot(diffuse_median,corrected_tilted_diffuse, '.-', [0,500],[0,500],'r--');
    xlabel('diffuse median');
    ylabel('tilted diffuse');

    figno=11;
    figure(figno); 
    ax(1) = subplot(2,4,3:4); % plot(times,ClearSkyGHI);
%     plot(times,irr,'r.-');
    plot(times,100*(diffuse_median-tilted_diffuse)./diffuse_median,'g.-');
    grid on;
    datetickzoom;
    title('total irradiation');
    ax(2) = subplot(2,4,7:8); plot(times,diffuse_median,'b.-',times,tilted_diffuse,'r.--');
    grid on;
    datetickzoom;
    title('diffuse irradiation');
    linkaxes(ax,'x');
    show_upd(1);
    datacursormode on;
    dcmobj = datacursormode(gcf); 
    set(dcmobj,'UpdateFcn',@show_upd,'DisplayStyle','window');
    function out = show_upd(i, event_obj)
        if nargin==2, i=event_obj.DataIndex; end
        d = data{i};
        figure(figno);
        subplot(2,4,[1 2 5 6]); 
        imshow([img_save_path d.day '__' num2str(d.j) '.jpeg']);
        drawnow;
        if nargin<2, out=[]; 
        else
            out = {datestr(event_obj.Position(1),'yyyy-mm-dd HH:MM:SS'), ['clouds: ',num2str(d.clouds),'%'], ['cloud-shine: ', num2str(d.cloud_shine)],['sun_flag: ',num2str(d.sun_flag)],['irr: ', num2str(d.irr(1))],['diffuse: ', num2str(d.diff_irr)], ['tilted_diffuse: ', num2str(tilted_diffuse(i))], ['clear_DNI: ', num2str(d.clear_irr(2))],['zenith: ',num2str(d.zenith)]};
        end
    end
end
