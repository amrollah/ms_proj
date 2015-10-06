function LabelerGUI_v2(labeler_struct_v2,varargin)

figure('KeyPressFcn',@cb)


%im_list = {'sky_image0.jpg','sky_image1.jpg','sky_image2.jpg','sky_image3.jpg','sky_image4.jpg','sky_image5.jpg'};
%im_arr = cellfun(@(x) sph2rect(double(imread(x))),im_list(:),'UniformOutput',false);
highResImage = [];
currentLabel = [];
lowResRect = [];
im_labelers = initializeLabelers();
clr_map = getColorMap();
lastPosition = [10 10 40 40];

im_str = labeler_struct_v2;
obj = setupSystemObject(im_str.im_dir);

nargs = numel(varargin);
if(numel(varargin)>0)
    im_labelers = varargin{1};
    obj.frameCount = length(im_labelers);
end

a1 = subplot(1,2,1);
a2 = subplot(1,2,2);
displayImage();


    function obj = setupSystemObject(image_folder)

        obj.reader = image_folder;
        obj.vid_ims = ls(obj.reader); % get the list of files in the folder
        obj.im_list = cellstr(obj.vid_ims(3:end,:)); % get the name list of those list of files
        obj.finalFrame = length(obj.im_list);
        
        obj.frameCount = 1; % the order of frames that are processed.
        obj.frameRate = 1;
        obj.frameName = []; % the file name of the frame that is processed the last.
    end

    function im_labelers = initializeLabelers()
        % create an empty array of tracks
        im_labelers = struct(...
            'id', {}, ...
            'name', {}, ...
            'im_labels', {},...
            'im_rects', {});
        
    end

    function im_labels = initializeLabels()
        % create an empty array of tracks
        im_labels = struct(...
            'id', {}, ...
            'name', {}, ...
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
                'id', obj.frameCount, ...
                'name', obj.frameName, ...
                'bbox', lowResRect{i}.getPosition(), ...
                'mask', lowResRect{i}.createMask(),...
                'color', lowResRect{i}.getColor()*255, ...
                'class', clr_map((lowResRect{i}.getColor()*255)*[1e6;1e3;1]));
            im_labels(i) = new_label;
            lastPosition = new_label.bbox;
        end
        new_labeler = struct(...
            'id', obj.frameCount, ...
            'name', obj.frameName,...
            'im_labels', {im_labels},...
            'im_rects', {lowResRect});
        
        im_labelers(obj.frameCount) = new_labeler;
    end

    function frame = readFrame()
        % in this function we read the frame on the order from the folder
        
        % if the destination buffer do not exist, than we immediately
        % terminate
        if(exist(obj.reader) ~=7)
            disp('Input stream should have been closed... Terminating...');
            frame = []; % indicator of an unsuccessful read trial
            return;
        end
        
        % if it is not a live operation, directly read the file in
        % order
        if(length(obj.im_list) >= obj.frameCount) % if there are still some images to be read
            obj.frameName = char(strcat(obj.reader, '/', obj.im_list(obj.frameCount)));
            frame = imread(obj.frameName);  % read the image directly, no need for the lock operations
        else
            frame = []; % if no more frames in the recording, directly terminate
            display('Done reading strored images... Terminating...');
        end
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
                    obj.frameCount = mod(obj.frameCount,obj.finalFrame)+obj.frameRate;
                    displayImage();
                case 'leftarrow'
                    recordLabels();
                    obj.frameCount = mod(obj.frameCount+obj.finalFrame-(2*obj.frameRate),obj.finalFrame)+obj.frameRate;
                    displayImage();
                case 'downarrow'
                    recordLabels();
                case 'b'
                    %current_folder = pwd;
                    %cd(im_str.rec_folder);
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
                    %cd(current_folder);
                    cd ../..
                otherwise
            end
        end
    end

    function displayImage()
        def_im = readFrame();
        highResImage = def_im;
        lowResImage = imresize(highResImage,1);
        imshow(lowResImage,'Parent',a1);
        addLabel(1);
    end

    function addLabel(init)
        if(init)
            currentLabel = 0;
            lowResRect = [];
            if(length(im_labelers)>=obj.frameCount)
                im_labels = im_labelers(obj.frameCount).im_labels;
                for i = 1:length(im_labels)
                    currentLabel = currentLabel + 1;
                    lowResRect{currentLabel} = imrect(a1,im_labels(currentLabel).bbox);
                    lastPosition = im_labels(currentLabel).bbox;
                    lowResRect{currentLabel}.setColor(im_labels(currentLabel).color/255);
                    lowResRect{currentLabel}.addNewPositionCallback(@(pos)Callback(pos,a2,highResImage));
                    fcn = makeConstrainToRectFcn('imrect',get(a1,'XLim'),get(a1,'YLim'));
                    setPositionConstraintFcn(lowResRect{currentLabel},fcn);
                    Callback(lowResRect{currentLabel}.getPosition(),a2,highResImage);
                end
                if(length(im_labels)~=0)
                    return
                end
            end
        end
        currentLabel = currentLabel + 1;
        initialPosition = lastPosition;
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