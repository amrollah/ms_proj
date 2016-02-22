function diffuse()
    img_save_path='';
    prj_path='';
    proj_path;
    load('calc\data_clean.mat', 'data');
    times = cellfun(@(d) d.time, data);
    irr = cellfun(@(d) d.irr(1), data);
    diffuse = cellfun(@(d) d.diff_irr, data);
    figno=11;
    figure(figno); 
    ax(1) = subplot(2,4,3:4); plot(times,irr,'r.-');
    grid on;
    datetickzoom;
    title('total irradiation');
    ax(2) = subplot(2,4,7:8); plot(times,diffuse,'b.-');
    grid on;
    datetickzoom;
    title('diffuse irradiation');
    linkaxes(ax,'x');
    show_upd(1);
    datacursormode on;
    dcmobj = datacursormode(gcf); 
    set(dcmobj,'UpdateFcn',@show_upd,'DisplayStyle','window');
    function out = show_upd(i, event_obj)
        if nargin<2, d = data{i};
        else d = data{event_obj.DataIndex}; end
        figure(figno);
        subplot(2,4,[1 2 5 6]); 
        imshow([img_save_path d.day '__' num2str(d.j) '.jpeg']);
        drawnow;
        if nargin<2, out=[]; 
        else
            out = {datestr(event_obj.Position(1),'yyyy-mm-dd HH:MM:SS'), ['clouds: ',num2str(d.clouds),'%'], ['cloud-shine: ', num2str(d.cloud_shine)],['sun_flag: ',num2str(d.sun_flag)],['irr: ', num2str(d.irr(1))],['diff_irr: ', num2str(d.diff_irr)],['clear_DNI: ', num2str(d.clear_irr(2))],['zenith: ',num2str(d.zenith)]};
        end
    end
end
