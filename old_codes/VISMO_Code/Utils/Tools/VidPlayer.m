%% This tool is prepared to be able to play videos in real time on matlab. 
% in a sense, this is a custom video player. 

% B. Zeydan, 27. Feb. 2013

classdef VidPlayer
%    UNTITLED Summary of this class goes here
%    Detailed explanation goes here
    
    properties (GetAccess=private)
        fig_handler
        ax_handler
    end
    
    methods
        function obj = VidPlayer(position,varargin)
            if(~isempty(varargin))
                buttonPressFunction = varargin(1);
                obj.fig_handler = figure('KeyPressFcn',buttonPressFunction);
            else
                obj.fig_handler = figure;
            end
            %box off; % for the initial frame, not to see the box, does not have any effect afterwards
            axis off; % for the initial frame, not to see the axis.  
            obj.ax_handler = axes;
            set(obj.ax_handler, 'visible', 'off')
            set(obj.fig_handler,'Position',position);
            set(obj.fig_handler,'Units','pixels');
        end
        
        function Step(obj,im,varargin)
            if(~isempty(varargin))
                titleText = varargin(1);
            else
                titleText = ' ';
            end
            %cla(obj.ax_handler);
%             imshow(im,'Parent',obj.ax_handler);
            if(min(im(:))<0)
                im = double(im) - min(im(:));
            end
            if(max(im(:))>1)
                im = double(im)/double(max(im(:)));
            end
            imagesc(im,'Parent',obj.ax_handler); 
            colormap(obj.ax_handler,gray); %may cause issues.
            axis (obj.ax_handler, 'image');
            set(obj.ax_handler,'xtick',[],'ytick',[]) %disable this if you
%             want to see the axes corresponding to pixels
            title(obj.ax_handler,titleText);
            
            drawnow;
        end
        
        function rgbframe = GetFrame(obj)
            F = getframe(obj.fig_handler);
            rgbframe = F.cdata;
        end
    end
    
end