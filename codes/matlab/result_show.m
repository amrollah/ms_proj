function result_show(data,xt,yt)
    close all;
    img_save_path='';
    prj_path='';
    proj_path;
    figno=11;
    figure(figno); 
    subplot(2,4,[3,4,7,8]);
    scatter(xt,yt,'b','.');
    lsline;
    hold on;
    plot([0 400], [0 400],'r-');
    xlabel('Predit Irradiance W/M^2');
    ylabel('Measured Irradiance W/M^2');
    xlim([0 400]);
    ylim([0 400]);
    title('SV ordinal regression machine');
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
            out = {datestr(d.time,'yyyy-mm-dd HH:MM:SS'), ['i: ', num2str(i)], ['clouds: ',num2str(d.clouds),'%'], ['sat_fact: ', num2str(d.sat_fact)],['sun_flag: ',num2str(d.sun_flag)],['irr: ', num2str(d.irr(1))],['diffuse: ', num2str(d.diff_irr)], ['tilted_diffuse: ', num2str(d.corr_tilt_diff)], ['clear_diffuse: ', num2str(d.clear_irr(3))],['zenith: ',num2str(d.zenith)],['cc: ' num2str(d.cc_fact)]};
        end
    end
end
