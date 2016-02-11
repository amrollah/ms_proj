function diffuse()
    load('img_data3.mat');
    if nargin<2, method = 'showtraj'; end
    figno=11;
    figure(figno); title('Diffuse vs Global irradiation');
    ax(1) = subplot(2,4,3:4); obj.plotdiffuse;
    title(obj.data.day,'interpreter','none');
    ax(2) = subplot(2,4,7:8); obj.plotirr;
    linkaxes(ax,'x');
    j=1;
    show_upd(j);
    dcmobj = datacursormode(gcf); 
    set(dcmobj,'UpdateFcn',@show_upd);
    function out = show_upd( j, event_obj)
        if nargin>=2
          [~,j] = min(abs(event_obj.Position(1)-obj.data.ti));
        end
        figure(figno);
        subplot(2,4,[1 2 5 6]); 
        obj.(method)(j);
        if ~isfield(obj.calc, 'seg') || length(obj.calc.seg.clear_sun_flag)<j || obj.calc.seg.clear_sun_flag(j)==0
            obj.loadframe(j);
            subplot(2,4,3:4); obj.plotdiffuse;
        end
        drawnow;
        if nargin<2, out=[]; 
        else out = {datestr(event_obj.Position(1),'HH:MM:SS'), ['clouds: ',num2str(obj.calc.seg.clouds_fact(j)),'%'], ['cloud-shine: ', num2str(obj.calc.seg.cloud_shine_fact(j))]};
        end
    end
end