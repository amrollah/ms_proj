function diffuse()
    img_save_path = 'C:\data\diffuse\';
    load('calc\img_data.mat', 'data');
    irr=
    figno=11;
    figure(figno); title('Diffuse vs Global irradiation');
    ax(1) = subplot(2,4,3:4); %obj.plotdiffuse;
    title(obj.data.day,'interpreter','none');
    ax(2) = subplot(2,4,7:8); %obj.plotirr;
    linkaxes(ax,'x');
    i=1;
    show_upd(i);
    dcmobj = datacursormode(gcf); 
    set(dcmobj,'UpdateFcn',@show_upd);
    function out = show_upd( i, event_obj)
        d = data{i};
        figure(figno);
        subplot(2,4,[1 2 5 6]); 
        imshow([img_save_path d.day '__' num2str(d.j) '.jpeg']);
        drawnow;
        if nargin<2, out=[]; 
        else out = {datestr(event_obj.Position(1),'yyyy-mmmm-dd HH:MM:SS'), ['clouds: ',num2str(d.clouds_fact),'%'], ['cloud-shine: ', num2str(d.cloud_shine_fact)]};
        end
    end
end