function LabelerGUI(labeler_struct,varargin)

figure('KeyPressFcn',@cb)


%im_list = {'sky_image0.jpg','sky_image1.jpg','sky_image2.jpg','sky_image3.jpg','sky_image4.jpg','sky_image5.jpg'};
%im_arr = cellfun(@(x) sph2rect(double(imread(x))),im_list(:),'UniformOutput',false);
currentFrame = 1;
highResImage = [];
currentLabel = [];
lowResRect = [];
im_labelers = initializeLabelers();
clr_map = getColorMap();

im_arr = labeler_struct.im_arr;
im_str = labeler_struct.im_str;

nargs = numel(varargin);
if(numel(varargin)>0)
    im_labelers = varargin{1};
    currentFrame = length(im_labelers);
end

a1 = subplot(1,2,1);
a2 = subplot(1,2,2);
displayImage();


    function im_labelers = initializeLabelers()
        % create an empty array of tracks
        im_labelers = struct(...
            'id', {}, ...
            'im_labels', {},...
            'im_rects', {});
        
    end

    function im_labels = initializeLabels()
        % create an empty array of tracks
        im_labels = struct(...
            'id', {}, ...
            'bbox', {}, ...
            'mask', {}, ...
            'color', {}, ...
            'class', {});
        
    end

    function recordLabels()
        % create an empty array of labels
        im_labels = initializeLabels();
        
        for i = 1:length(lowResRect)
            new_label = struct(...
                'id', currentFrame, ...
                'bbox', lowResRect{i}.getPosition(), ...
                'mask', lowResRect{i}.createMask(),...
                'color', lowResRect{i}.getColor()*255, ...
                'class', clr_map((lowResRect{i}.getColor()*255)*[1e6;1e3;1]));
            im_labels(i) = new_label;
        end
        new_labeler = struct(...
            'id', currentFrame, ...
            'im_labels', {im_labels},...
            'im_rects', {lowResRect});
        
        im_labelers(currentFrame) = new_labeler;
    end

    function Callback(position,axesHandle, highResImage)
        
        position = position * 1;
        x1 = position(1);
        y1 = position(2);
        x2 = position(1) + position(3);
        y2 = position(2) + position(4);
        
        %further position constraints added for all dimensions
        if(x1<1)
            x1 = 1;
        end
        if(x2 > size(highResImage,2))
            x2 = size(highResImage,2);
        end
        if(y1<1)
            y1 = 1;
        end
        if(y2 > size(highResImage,1))
            y2 = size(highResImage,1);
        end
        
        highResThumbnail = highResImage(round(y1:y2),round(x1:x2),:);
        
        if isempty( get(axesHandle,'Children'))
            imshow(highResThumbnail,'Parent',axesHandle);
        else
            imHandle = get(axesHandle,'Children');
            oldSize = size(get(imHandle,'CData'));
            if ~isequal(oldSize, size(highResThumbnail))
                imshow(highResThumbnail,'Parent',axesHandle);
            else
                set( imHandle,'CData', highResThumbnail);
            end
        end
    end

    function cb(~,evnt)
        if strcmp(evnt.Modifier{1},'control')
            switch evnt.Key
                case 'uparrow'
                    addLabel(0);
                case 'rightarrow'
                    recordLabels();
                    currentFrame = mod(currentFrame,length(im_arr))+1;
                    displayImage();
                case 'leftarrow'
                    recordLabels();
                    currentFrame = mod(currentFrame+length(im_arr)-2,length(im_arr))+1;
                    displayImage();
                case 'downarrow'
                    recordLabels();
                case 'b'
                    e_val = exist('recordings');
                    if ~(e_val == 7)
                        mkdir recordings;
                    end
                    cd recordings;
                    c_list = ls;
                    b_num = sort(str2num(c_list(3:end,:)));
                    
                    if isempty(b_num)
                        b_num = 0;
                    else
                        b_num = b_num(end);
                    end
                    
                    b_num = b_num + 1;
                    b_str = num2str(b_num);
                    mkdir(b_str);
                    cd(b_str);
                    
                    save('labels', 'im_labelers');
                    save('details', 'im_str');
                    save_folder = pwd;
                    save('saved','save_folder');
                    cd ../..
                otherwise
            end
        end
    end

    function displayImage()
        def_im = im_arr{currentFrame};
        highResImage = def_im;
        lowResImage = imresize(highResImage,1);
        imshow(lowResImage,'Parent',a1);
        addLabel(1);
    end

    function addLabel(init)
        if(init)
            currentLabel = 0;
            lowResRect = [];
            if(length(im_labelers)>=currentFrame)
                im_labels = im_labelers(currentFrame).im_labels;
                for i = 1:length(im_labels)
                    currentLabel = currentLabel + 1;
                    lowResRect{currentLabel} = imrect(a1,im_labels(currentLabel).bbox);
                    lowResRect{currentLabel}.setColor(im_labels(currentLabel).color/255);
                    lowResRect{currentLabel}.addNewPositionCallback(@(pos)Callback(pos,a2,highResImage));
                    fcn = makeConstrainToRectFcn('imrect',get(a1,'XLim'),get(a1,'YLim'));
                    setPositionConstraintFcn(lowResRect{currentLabel},fcn);
                    Callback(lowResRect{currentLabel}.getPosition(),a2,highResImage);
                end
                return
            end
        end
        currentLabel = currentLabel + 1;
        initialPosition = [10 10 40 40];
        lowResRect{currentLabel} = imrect(a1,initialPosition);
        lowResRect{currentLabel}.addNewPositionCallback(@(pos)Callback(pos,a2,highResImage));
        fcn = makeConstrainToRectFcn('imrect',get(a1,'XLim'),get(a1,'YLim'));
        setPositionConstraintFcn(lowResRect{currentLabel},fcn);
        Callback(initialPosition, a2, highResImage);
    end

    function clr_map = getColorMap()
        colors_rgb = [248 72 248; ...
            72 248 72; ...
            248 246 74; ...
            72 136 248; ...
            72 72 248; ...
            248 79 79; ...
            0 0 0; ...
            232 232 232; ...
            72 248 248];
        colors_keys = 1e6*colors_rgb(:,1)+1e3*colors_rgb(:,2)+colors_rgb(:,3);
        
        colors_class = {'m', 'g', 'y', 'b', 'p', 'r', 'k', 'w', 'c'};
        clr_map = containers.Map(colors_keys(:),colors_class(:));
    end
end