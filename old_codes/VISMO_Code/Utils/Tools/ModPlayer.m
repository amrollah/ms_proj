%% This tool is prepared to be able to play the model video in real time on matlab.
% in a sense, this is a custom model player.

% B. Zeydan, 12. Jun. 2013
% B. Zeydan, 03. Jul. 2013 , minor update for the plotting functionality
classdef ModPlayer
    
    properties (GetAccess=private)
        fig_handler
        ax_handler
        leg_handler
        object_handlers
        option
        coords
        image
        image_RealSize
        image_MapFrame
    end
    
    methods
        function obj = ModPlayer(position,option,varargin)
            
            obj.fig_handler = figure;
            obj.ax_handler = axes;
            obj.leg_handler = legend(obj.ax_handler);
            %axis off;
            %box off;
            obj.option = option;
            obj.object_handlers = [];
            if(~isempty(varargin))
                obj.image = varargin{1};
                if(length(varargin)>=2) %get the coords
                    obj.coords = varargin{2};
                    if(length(varargin)>=3) %get the realsize image
                        obj.image_RealSize = imread(varargin{3});
                    end
                end
            end
            switch option
                case 'rectangle'
                    set(obj.ax_handler,'XLim',[-170000,170000]);
                    set(obj.ax_handler,'YLim',[-170000,170000]);
                case 'circle'  
                    set(obj.ax_handler,'XLim',[-170000,170000]);
                    set(obj.ax_handler,'YLim',[-170000,170000]);
                    hold on;
                case 'plot'
                    %set(obj.ax_handler,'XLimMode','auto');
                    %set(obj.ax_handler,'YLimMode','auto');
                    %set(obj.ax_handler,'XLimMode','manual');
                    %set(obj.ax_handler,'YLimMode','manual');
                    set(obj.ax_handler,'XLim',[0,15]);
                    set(obj.ax_handler,'YLim',[0,1300]);
                    set(obj.ax_handler,'FontSize',13);
                case 'generic'
                    
                case 'map'
                    if(~isempty(varargin))
                        CHCRC = varargin{1};
                        coords = obj.coords;
                        x_lb = coords(1);
                        x_ub = coords(2);
                        y_lb = coords(3);
                        y_ub = coords(4);
                        if(length(coords)>6)
                            x_s = coords(7);
                            y_s = coords(8);
                        end
                        surface('XData',[x_lb x_ub; x_lb x_ub],...
                                    'YData',[y_lb y_lb; y_ub y_ub],...
                                    'ZData',[0 0; 0 0],'CData',flipdim(CHCRC,1),...
                                    'FaceColor','texturemap','EdgeColor','none');
                        hold on; 
                        plot([0 0],[y_lb y_ub],'m') % draw y axis
                        plot([x_lb x_ub],[0 0],'m') % draw x axis
                        
                        if(length(coords)>6)
                            scatter3(x_s,y_s,0,100,'r','fill');
                        end
                        hold off;
                        xlabel(obj.ax_handler,'--X axis--');
                        ylabel(obj.ax_handler,'--Y axis--');
                        zlabel(obj.ax_handler,'--Z axis--');
%                         view(-20,49);
                    end
                case 'projected map'
                    if(~isempty(varargin))
                        CHCRC = varargin{1};
                        coords = obj.coords;
                        x_lb = coords(1);
                        x_ub = coords(2);
                        y_lb = coords(3);
                        y_ub = coords(4);
                        if(length(coords)>6)
                            x_s = coords(7);
                            y_s = coords(8);
                        end
                        imagesc([x_lb x_ub], [y_lb y_ub], CHCRC, 'Parent',obj.ax_handler);
                        set(obj.ax_handler,'xtick',[],'ytick',[])
                        axis (obj.ax_handler, 'image');
                        xlabel(obj.ax_handler,'--X axis--');
                        ylabel(obj.ax_handler,'--Y axis--');
%                         zlabel(obj.ax_handler,'--Z axis--');
%                         view(0,90);
                    end
                otherwise
            end
            set(obj.fig_handler,'Position',position);
            set(obj.fig_handler,'Units','pixels');
        end
        
        function obj = Step(obj,boxes,varargin)
            if(~isempty(varargin))
                if(length(varargin)>3)
                    XLabelText = varargin(1);
                    YLabelText = varargin(2);
                    titleText  = varargin(3);
                    legendText = cellstr(varargin(4));
                elseif(length(varargin)>2)
                    XLabelText = varargin(1);
                    YLabelText = varargin(2);
                    titleText  = varargin(3);
                    legendText = [];
                end
            end
            
            for i = 1:length(obj.object_handlers)
                delete(obj.object_handlers(i));
            end
            obj.object_handlers = [];
            
            switch obj.option
                case 'reactangle'
                    for i = 1:size(boxes,1)
                        %rectangle('Position',[boxes(i,2),boxes(i,1),boxes(i,4),boxes(i,3)],'Parent',obj.ax_handler);
                        if(boxes(i,3)>0 && boxes(i,4)>0)
                            obj.object_handlers(i)=rectangle('Position',boxes(i,:),'Parent',obj.ax_handler);
                        end
                        %rectangle('Position',boxes(i,:));
                    end
                    drawnow;
                case 'circle'
                    
                    for i = 1:size(boxes,1)
                        if(boxes(i,3)>0)
                            % faster this way
                            x = boxes(i,1); y = boxes(i,2); r = boxes(i,3);
                            ang=0:0.01:2*pi;
                            xp=r*cos(ang);
                            yp=r*sin(ang);
                            obj.object_handlers(i)=plot(x+xp,y+yp,'Parent',obj.ax_handler);
                            %obj.circle(boxes(i,1),boxes(i,2),boxes(i,3));
                        end
                    end
                    drawnow;
                case 'plot'
                    %cla(obj.ax_handler,'reset');
                    %set(obj.ax_handler,'XLimMode','auto');
                    %set(obj.ax_handler,'YLimMode','auto');
                    
                    if(size(boxes,2)==2)
                        obj.object_handlers = plot(boxes(:,1),boxes(:,2),'Parent',obj.ax_handler,'LineWidth',4);
                    elseif(size(boxes,2)==4)
                        obj.object_handlers = plot(boxes(:,1),boxes(:,2),boxes(:,3),boxes(:,4),'Parent',obj.ax_handler,'LineWidth',4);
                    else
                        obj.object_handlers = plot(boxes(:,1),boxes(:,2),boxes(:,3),boxes(:,4),boxes(:,5),boxes(:,6),'Parent',obj.ax_handler,'LineWidth',4);
                    end
                    
                    set(obj.ax_handler,'XLim',[min(boxes(:,1)),max(boxes(:,1))]);
                    if(boxes(:,4)<5)
                        set(obj.ax_handler,'YLim',[0,1]);
                    else
                        set(obj.ax_handler,'YLim',[0,1300]);
                    end
                    
                    xlabel(obj.ax_handler,XLabelText);
                    ylabel(obj.ax_handler,YLabelText);
                    title(obj.ax_handler,titleText);
                    if(~isempty(legendText))
                        legend(obj.object_handlers,legendText,'Location','Southeast');
                    end
                    drawnow;
                case 'generic'
                    cla(obj.ax_handler,'reset');
                    plot(boxes(:,1),boxes(:,2),'Parent',obj.ax_handler,'LineWidth',4);
                    xlabel(obj.ax_handler,XLabelText);
                    ylabel(obj.ax_handler,YLabelText);
                    title(obj.ax_handler,titleText);
                    if(~isempty(legendText))
                        legend(obj.object_handlers,legendText,'Location','best');
                    end
                    drawnow;
                case 'map'
                    if(isempty(varargin))
                        colors = 0.2*ones(size(boxes));
                    else
                        if(sum(size(varargin{1})-size(boxes))~=0)
                            disp('every patch has to have a color...');
                        else
                            colors = varargin{1};
                        end
                    end
                            
                    for i = 1:size(boxes,2)
                        contourPts = vertcat(boxes{1,i}{1,:});
                        x_cloud = contourPts(:,1)';
                        y_cloud = contourPts(:,2)';
                        z_cloud = contourPts(:,3)';
                        obj.object_handlers(i) = patch ( x_cloud, y_cloud, z_cloud, colors(i),'Parent',obj.ax_handler);
                        alpha (obj.object_handlers(i),0.9);
                    end
                case 'projected map'
                    if(isempty(varargin))
                        colors = 0.2*ones(size(boxes));
                    else
%                         if(sum(size(varargin{1})-size(boxes))~=0)
%                             disp('every patch has to have a color...');
%                         else
%                             colors = varargin{1};
%                         end
                    end
                    coords = obj.coords;
                    x_lb = coords(1);
                    x_ub = coords(2);
                    y_lb = coords(3);
                    y_ub = coords(4);
                    im_dims = size(obj.image);
                    tempImg = obj.image;
                    for contourInd = 1:length(boxes)
                        for i = 1:size(boxes{contourInd},2)
                            contourPts = vertcat(boxes{contourInd}{1,i}{1,:});
                            x_cloud = contourPts(:,1);
                            y_cloud = contourPts(:,2);
                            %%resize the contours so that they fit the resized
                            %%image and are transformed correctly between the two
                            %%mappings of (image) and (real world)
                            x_cloud = ((x_cloud-x_lb)./(x_ub-x_lb).*(im_dims(2)))+1;
                            y_cloud = ((-(y_cloud-y_lb)./(y_ub-y_lb))+1).*(im_dims(1))+1;
                            
                            
                            
                            predShadowContoursOnMap{1,i} = num2cell([x_cloud, y_cloud],2)';
                            %                         obj.object_handlers(i) = patch ( x_cloud, y_cloud, z_cloud, colors(i),'Parent',obj.ax_handler);
                            %                         alpha (obj.object_handlers(i),0.9);
                            
                        end
                        emphasisOnClouds = 0.5; %should be between [0,1]; 1 being opaque.
                        if (length(varargin{1}{contourInd})==3) % if it's the color info
                            colors = varargin{1}{contourInd};
                            mapProj = cv.drawContours(tempImg, predShadowContoursOnMap,'Color', colors, 'Thickness', -1);
                            mapProj = uint8(double(mapProj)*emphasisOnClouds+ double(tempImg)*(1-emphasisOnClouds));
                        else
                            mapProj = cv.drawContours(tempImg, predShadowContoursOnMap,'Color', [255 0 0], 'Thickness', -1);
                            mapProj = uint8(double(mapProj)/2+ double(tempImg)/2);
                        end
                        tempImg = mapProj;
                    end
                    imagesc([x_lb x_ub], [y_lb y_ub],mapProj, 'Parent',obj.ax_handler); % Add the axis info here from obj or model3D
                    title(obj.ax_handler,'Experimental ID');
                    set(obj.ax_handler,'xtick',[],'ytick',[]);
                    axis (obj.ax_handler, 'image', 'tight');
%                     drawnow;
                    obj.image_MapFrame=mapProj;
                    imwrite(mapProj, [char(varargin{2}) '_map.jpeg']); % delete this later.
                otherwise
                    disp('please enter a valid option such as circle or rectangle...');
            end
        end
%         function rgbframe = getMapFrame(obj) %%FIX: image_MapFrame
%         returns empty.
%             rgbframe = obj.image_MapFrame;
%             if(isempty(rgbframe))
%                 rgbframe = obj.image;
%             end
%             
%         end
        
        function rgbframe = GetFrame(obj)
%             movegui(obj.fig_handler);
%             F = getframe(obj.fig_handler);
%             rgbframe = F.cdata;
              rgbframe = obj.image;
        end
        
        function circle(obj,x,y,r)
            %x and y are the coordinates of the center of the circle
            %r is the radius of the circle
            %0.01 is the angle step, bigger values will draw the circle faster but
            %you might notice imperfections (not very smooth)
            ang=0:0.01:2*pi;
            xp=r*cos(ang);
            yp=r*sin(ang);
            plot(x+xp,y+yp,'Parent',obj.ax_handler);
        end
    end
    
end