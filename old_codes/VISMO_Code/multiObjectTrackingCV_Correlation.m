%% Cloud tracking using correlation method for speed estimation.
%
% Multiple object tracking with correlation for could tracking:
%
% video_name: path to the input video file
% record: option indicating that a record should be taken or not 
% verbose: option indication that the video should be seen or not
% finalFrame: number of frames of the video, that the algorithm will run
% on.
% B. Zeydan, 27. Feb. 2013


function [recording_obs,recording_pred] = multiObjectTrackingCV_Correlation(image_folder,record,verbose,finalFrame,imageMask,sunPattern,calibModel,ecalibModel)
% create system objects used for reading video, detecting moving objects,
% and displaying the results

if (nargin < 2)
    record = 0;
end

if(nargin < 3)
    verbose = 1;
end

if(nargin < 4)
    finalFrame = 10000;
end

%video_name = 'sky_low.mp4'; % the name of the stream file (video or a collection of images)
obj = setupSystemObjects(image_folder,verbose,imageMask,calibModel,ecalibModel);
scale_factor = 0.25;
sunPattern = imresize(sunPattern,scale_factor);
ff = readFrame();
obj.frameCount = obj.frameCount-obj.frameRate;
frameSize = size(ff);
%frameSize = [frameSize(2),frameSize(1)];
R = 10;
y_edges = [1:R:frameSize(1)-R+1];
x_edges = [1:R:frameSize(2)-R+1];
[bby,bbx] = meshgrid(y_edges,x_edges);

mesh_lut_temp = [bbx(:),bby(:),(R-1)*ones(size([bbx(:),bby(:)]))];
mesh_cent_lut_temp = [bbx(:),bby(:)] + 0.5*R*ones(size([bbx(:),bby(:)]));
iter2 = 1;
mesh_lut = [];
mesh_cent_lut = [];

for iter = 1:size(mesh_cent_lut_temp)
    if(obj.imageMask(round(mesh_cent_lut_temp(iter,2)./scale_factor),round(mesh_cent_lut_temp(iter,1)./scale_factor),1))
        mesh_lut(iter2,:) = mesh_lut_temp(iter,:);
        mesh_cent_lut(iter2,:) = mesh_cent_lut_temp(iter,:);
        iter2 = iter2+1;
    end
end

sunEliminationCycle = 100; % every this many frames, we calculate the sun's position once again. 
sunPosition = [-1,-1]; % we initialize the position of sun, to be used throughout the code

tracks = initializeTracks(); % create an empty array of tracks
nextId = 1; % ID of the next track
areas = [];
centroids = [];
bboxes = [];
mask = [];

maskMethod = 'body';
trackOption = 'speed';

if(record && (predictionCycle > updateCycle))
    tracks2 = initializeTracks(); % create an empty array of tracks
    nextId2 = 1; % ID of the next track
    trackOption2 = 'pose_v';
    areas2 = [];
    centroids2 = [];
    bboxes2 = [];
    mask2 = [];
end

if(record)
    recording_pred = cell(finalFrame,1);
    recording_obs = cell(finalFrame,1);
end

trackOption = 'pose_v';
areas = [];
centroids = [];
bboxes = [];
mask = [];
cloud_mask = [];
frame = [];

% detect moving objects, and track them across video frames
while obj.frameCount < obj.finalFrame %~isDone(obj.reader)
    tic
    prev_frame = frame;
    frame = readFrame();
    if((sum(sunPosition)<0) || mod(obj.frameCount,sunEliminationCycle)>=(sunEliminationCycle-obj.frameRate))
        sunPosition = getSunPosition(frame,sunPattern);
        sunPatch = getSunPatch(sunPosition);
        sunPatch = imresize(sunPatch,scale_factor);
    end
    
    prev_mask = mask;
    prev_cloud_mask = cloud_mask;
    
    frame = eliminateSun(frame,sunPosition);
    %[mask,cloud_mask] = getBinaryImage(frame,maskMethod,sunPatch);
    [mask,cloud_mask] = getBinaryImage2(frame,maskMethod);
    %mask = medfilt2(mask,[15 15]);

    if(obj.frameCount <= obj.frameRate+1)
        prev_frame = frame;
        prev_mask = mask;
        prev_cloud_mask = cloud_mask;
    end
    
    %tolmask = zeros(size(mask)+6*[R,R]);
    %tolmask(3*R:3*R+size(mask,1),3*R:3*R+size(mask,2)) = mask;
    %[max_val, max_ind] = cellfun(@(x) max(abs(normxcorr2(prev_mask(x(2):x(2)+x(4),x(1):x(1)+x(3))+mild_noise,tolmask(3*R+x(2)-3*R:3*R+x(2)+x(4)+3*R,3*R+x(1)-3*R:3*R+x(1)+x(3)+3*R)))) ,mat2cell(mesh_lut,ones(size(mesh_lut,1),1)), 'UniformOutput',false);
    %[max_val, max_ind] = cellfun(@(x) max(abs(normxcorr2(prev_mask(x(2):x(2)+x(4),x(1):x(1)+x(3)),mask))) ,mat2cell(mesh_lut,ones(size(mesh_lut,1),1)), 'UniformOutput',false);
    %[max_ind_x,max_ind_y] = cellfun(@(x) ind2sub(frameSize,x(1)),max_ind,'UniformOutput',false);
    
    [max_ind_x,max_ind_y] = cellfun(@(x) maxCorrPatchBackward(mask,cloud_mask,prev_mask,x),mat2cell(mesh_lut,ones(size(mesh_lut,1),1)),'UniformOutput',false);
    ll = [max_ind_y{:};max_ind_x{:}];
    %arrows = [mesh_cent_lut(:,2),mesh_cent_lut(:,1),ll(2,:)',ll(1,:)']-1;
    arrows = [mesh_cent_lut,ll'];
    predictionImage = uint8(zeros(size(frame)));
    cellfun(@(x) generatePrediction(frame,x),mat2cell([mesh_lut,ll'],ones(size([mesh_lut,ll'],1),1)),'UniformOutput',false);
    
    %flow = cv.calcOpticalFlowFarneback(prev_mask,mask);
    %[lx,ly] = meshgrid([1:size(flow,1)],[1:size(flow,2)]);
    %f1 = flow(:,:,1);
    %f2 = flow(:,:,2);
    %arrows = [lx(:),ly(:),f1(:),f2(:)];
    %if(~isempty(flow))
    %    figure, imshow(flow);
    %end
    if(verbose)
        displayTrackingResults();
    end
    toc
end

    function [max_ind_x, max_ind_y] = maxCorrPatchForward(patternImage,cloudMaskedPatternImage,searchImage,patchArea)
        searchArea = [patchArea(1)-0.3*patchArea(3),patchArea(2)-0.3*patchArea(4),patchArea(1)+1.3*patchArea(3),patchArea(2)+1.3*patchArea(4)];
        
        if(searchArea(1)<1)
            searchArea(1)=1;
        end
        if(searchArea(2)<1)
            searchArea(2)=1;
        end
        if(searchArea(3)>size(searchImage,2))
            searchArea(3)=size(searchImage,2);
        end
        if(searchArea(4)>size(searchImage,1))
            searchArea(4)=size(searchImage,1);
        end
        
        
        x = round(patchArea);
        y = round(searchArea);
        
        mild_noise = 0.0001*rand(R,R);
        
        patternIm = patternImage(x(2):x(2)+x(4),x(1):x(1)+x(3));
        patternCloudIm = cloudMaskedPatternImage(x(2):x(2)+x(4),x(1):x(1)+x(3));
        searchImPatternArea = searchImage(x(2):x(2)+x(4),x(1):x(1)+x(3));
        searchIm = searchImage(y(2):y(4),y(1):y(3));
        meanPatternIm = mean(patternIm(:));
        meanSearchImPatternArea = mean(searchImPatternArea(:));
        meanPatternCloudiness = mean(patternCloudIm(:));
        
        if((abs(meanPatternIm-meanSearchImPatternArea)>0.1*meanPatternIm)&&(meanPatternCloudiness>0.3))
            cc = normxcorr2(patternIm+mild_noise,searchIm);
            [max_val,max_ind] = max(abs(cc(:)));
            [max_ind_x,max_ind_y] = ind2sub(size(cc),max_ind(1));
            max_ind_x = max_ind_x-0.5*x(4)+y(2);
            max_ind_y = max_ind_y-0.5*x(3)+y(1);
            %max_ind_x = max_ind_x + y(2);
            %max_ind_y = max_ind_y + y(1);
        else
            max_ind_x = x(2)+(0.5*x(4));
            max_ind_y = x(1)+(0.5*x(3));
        end
    end

    function [max_ind_x, max_ind_y] = maxCorrPatchBackward(patternImage,cloudMaskedPatternImage,searchImage,patchArea)
        searchArea = [patchArea(1)-0.4*patchArea(3),patchArea(2)-0.4*patchArea(4),patchArea(1)+1.4*patchArea(3),patchArea(2)+1.4*patchArea(4)];
        
        if(searchArea(1)<1)
            searchArea(1)=1;
        end
        if(searchArea(2)<1)
            searchArea(2)=1;
        end
        if(searchArea(3)>size(searchImage,2))
            searchArea(3)=size(searchImage,2);
        end
        if(searchArea(4)>size(searchImage,1))
            searchArea(4)=size(searchImage,1);
        end
        
        x = round(patchArea);
        y = round(searchArea);
        
        mild_noise = 0.0001*rand(R,R);
        
        patternIm = patternImage(x(2):x(2)+x(4),x(1):x(1)+x(3));
        patternCloudIm = cloudMaskedPatternImage(x(2):x(2)+x(4),x(1):x(1)+x(3));
        searchImPatternArea = searchImage(x(2):x(2)+x(4),x(1):x(1)+x(3));
        searchIm = searchImage(y(2):y(4),y(1):y(3));
        meanPatternIm = mean(patternIm(:));
        meanSearchImPatternArea = mean(searchImPatternArea(:));
        meanPatternCloudiness = mean(patternCloudIm(:));
        
        if((abs(meanPatternIm-meanSearchImPatternArea)>0.1*meanPatternIm)||(meanPatternCloudiness>0.1))
            cc = normxcorr2(patternIm+mild_noise,searchIm);
            [max_val,max_ind] = max(cc(:));
            [max_ind_x,max_ind_y] = ind2sub(size(cc),max_ind(1));
            max_ind_x = max_ind_x-0.5*x(4)+y(2);
            max_ind_y = max_ind_y-0.5*x(3)+y(1);
            max_ind_x = 2*(x(2)+0.5*x(4)) - max_ind_x;
            max_ind_y = 2*(x(1)+0.5*x(3)) - max_ind_y;
            %max_ind_x = max_ind_x + y(2);
            %max_ind_y = max_ind_y + y(1);
        else
            max_ind_x = x(2)+(0.5*x(4));
            max_ind_y = x(1)+(0.5*x(3));
        end
    end
    

    function [max_ind_x, max_ind_y] = maxCorrPatchBackward2(patternImage,cloudMaskedPatternImage,searchImage,patchArea)
        searchArea = [patchArea(1)-0.4*patchArea(3),patchArea(2)-0.4*patchArea(4),patchArea(1)+1.4*patchArea(3),patchArea(2)+1.4*patchArea(4)];
        
        if(searchArea(1)<1)
            searchArea(1)=1;
        end
        if(searchArea(2)<1)
            searchArea(2)=1;
        end
        if(searchArea(3)>size(searchImage,1))
            searchArea(3)=size(searchImage,1);
        end
        if(searchArea(4)>size(searchImage,2))
            searchArea(4)=size(searchImage,2);
        end
        
        x = round(patchArea);
        y = round(searchArea);
        
        mild_noise = 0.0001*rand(R,R);
        
        patternIm = patternImage(x(1):x(1)+x(3),x(2):x(2)+x(4));
        patternCloudIm = cloudMaskedPatternImage(x(1):x(1)+x(3),x(2):x(2)+x(4));
        searchImPatternArea = searchImage(x(1):x(1)+x(3),x(2):x(2)+x(4));
        searchIm = searchImage(y(1):y(3),y(2):y(4));
        meanPatternIm = mean(patternIm(:));
        meanSearchImPatternArea = mean(searchImPatternArea(:));
        meanPatternCloudiness = mean(patternCloudIm(:));
        
        if((abs(meanPatternIm-meanSearchImPatternArea)>0.1*meanPatternIm)||(meanPatternCloudiness>0.1))
            cc = normxcorr2(patternIm+mild_noise,searchIm);
            [max_val,max_ind] = max(cc(:));
            [max_ind_x,max_ind_y] = ind2sub(size(cc),max_ind(1));
            max_ind_x = max_ind_x-0.5*x(3)+y(1);
            max_ind_y = max_ind_y-0.5*x(4)+y(2);
            max_ind_x = 2*(x(1)+0.5*x(3)) - max_ind_x;
            max_ind_y = 2*(x(2)+0.5*x(4)) - max_ind_y;
            %max_ind_x = max_ind_x + y(2);
            %max_ind_y = max_ind_y + y(1);
        else
            max_ind_x = x(1)+(0.5*x(3));
            max_ind_y = x(2)+(0.5*x(4));
        end
    end

    function generatePrediction(patternImage,patchArea_displacement)
        patchArea = patchArea_displacement(1:4);
        displacement = patchArea_displacement(5:6)-patchArea_displacement(3:4)/2;
        resultArea = [displacement(1),displacement(2),displacement(1)+patchArea(3),displacement(2)+patchArea(4)];
        %resultArea = [patchArea(1)+displacement(1),patchArea(2)+displacement(2),patchArea(1)+displacement(1)+patchArea(3),patchArea(2)+displacement(2)+patchArea(4)];
        
        if(resultArea(1)<1)
            patchArea(1) = patchArea(1)-resultArea(1);
            resultArea(1)=1;
            
        end
        if(resultArea(2)<1)
            patchArea(2) = patchArea(2)-resultArea(2);
            resultArea(2)=1;
        end
        if(resultArea(3)>size(predictionImage,2))
            patchArea(3) = patchArea(3)+size(predictionImage,2)-resultArea(3);
            resultArea(3) = size(predictionImage,2);
        end
        if(resultArea(4)>size(predictionImage,1))
            patchArea(4) = patchArea(4)+size(predictionImage,1)-resultArea(4);
            resultArea(4)=size(predictionImage,1);
        end
        
        x = round(patchArea);
        y = round(resultArea);
        
        mild_noise = 0.0001*rand(R,R);
        
        patternIm = patternImage(x(2):x(2)+x(4),x(1):x(1)+x(3),:);
        predictionImage(y(2):y(4),y(1):y(3),:) = predictionImage(y(2):y(4),y(1):y(3),:) + patternIm;
    end

    function obj = setupSystemObjects(image_folder,verbose,imageMask,calibModel,ecalibModel)
        % Initialize Video I/O
        % Create objects for reading a video from a file, drawing the tracked
        % objects in each frame, and playing the video.
        
        obj.live = 0;
        % if the read functionality is live, than we first check if the
        % buffer folder exists so the camera is already operational
        % if it cannot find the buffer folder, then it runs a bussy-wait
        % procedure on it. When it times out, it terminates the setup
        % function and terminates the program before ever proceeding. if
        % the read functionality is not live, then is is directly
        % initiated to be checked during the first read.
        if(obj.live)
            while(1)
                try_ct = try_ct + 1;
                if(exist(char([image_folder, '/Buffer']))==7)
                    obj.reader = char([image_folder, '/Buffer']);
                    break;
                else
                    if(try_ct > 10)
                        obj.reader = [];
                        disp('No input stream found.. Exit.');
                        return;
                    else
                        disp(char(['Trying to read for: ' num2str(try_ct)]));
                        pause(2);
                    end
                end
            end
            obj.vid_ims = [];
            obj.im_list = [];
            obj.finalFrame = inf;
        else
            obj.reader = image_folder;
            obj.vid_ims = ls(obj.reader); % get the list of files in the folder
            obj.im_list = cellstr(obj.vid_ims(3:end,:)); % get the name list of those list of files
            obj.finalFrame = length(obj.im_list);
        end
        obj.frameCount = 1; % the order of frames that are processed.
        obj.frameName = []; % the file name of the frame that is processed the last.
        obj.imageMask = imageMask;
        obj.frameRate = 2;
        obj.predictionImage = [];
        obj.calibModel = calibModel; % the internal calibration model that is used by the algorithm
        obj.ecalibModel = ecalibModel; % the external calibration model that is used by the algorithm
        imageMask = ones(obj.calibModel.height,obj.calibModel.width);
        [X,Y] = meshgrid(1:obj.calibModel.height,1:obj.calibModel.width);
        obj.X = X;
        obj.Y = Y;
        P = [X(:)';Y(:)'];
        p = cam2world(P,obj.calibModel);
        p = obj.ecalibModel.R*p;
        maskTemp = p(3,:)<-0.15;
        indicesTemp = sub2ind([obj.calibModel.height,obj.calibModel.width],P(1,:),P(2,:));
        imageMask(indicesTemp) = maskTemp;
        imageMask = logical(repmat(imageMask,[1,1,3]));
        obj.imageMask = imageMask;
        % create two video players, one to display the video,
        % and one to display the foreground mask
        if(verbose)
            obj.videoPlayer = VidPlayer([20, 400, 700, 400]);
            obj.maskPlayer = VidPlayer([740, 400, 700, 400]);
            obj.predPlayer = VidPlayer([740, 100, 700, 400]); % for the model in real dimensions
        end
    end

    function tracks = initializeTracks()
        % create an empty array of tracks
        tracks = struct(...
            'id', {}, ...
            'bbox', {}, ...
            'kalmanFilter', {}, ...
            'age', {}, ...
            'totalVisibleCount', {}, ...
            'consecutiveInvisibleCount', {},...
            'xspeed', {}, ...
            'yspeed', {}, ...
            'predSpeed', {}, ...
            'obsSpeed', {}, ...
            'area', {}, ...
            'speedFilter', {}, ...
            'positionFilter', {}, ...
            'centroidFilter', {}, ...
            'massFilter', {}, ...
            'centroid', {}, ...
            'predCentroid', {}, ...
            'prevCentroid', {});
        
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
        
        read_try = 1; % try to read a file for the first time
        if(obj.live) % if this is the live mode, meaning, we are going at the same paste as the camera
            vid_ims = ls(obj.reader); % get the list of files in the folder
            im_list = cellstr(vid_ims(3:end,:)); % get the name list of those list of files
            while(length(im_list{1}) == 0) % wait for an image to be written to the buffer
                if(exist(obj.reader) ~=7) % if the buffer folder is deleted during the wait, this means the camera has been shut down, so terminate
                    disp('Input stream should have been closed... Terminating...');
                    frame = []; % indicator of an unsuccessful read trial
                    return;
                end
                vid_ims = ls(obj.reader); % get the list of files in the folder
                im_list = cellstr(vid_ims(3:end,:)); % get the name list of those list of files
                read_try = read_try+1; % increment the number of trials and go for another try
                disp(char(['Trying to read for: ' num2str(read_try)])); % report the action to the user.
                if(read_try > 8) % if a read is tried for more than 8 times, then terminate with an unsuccessfull read
                    disp('Input stream not available... Terminating...');
                    frame = []; % unsuccessful read indicator
                    return;
                end
                pause(1); % wait 1 msec between two read trials, we can increase this time unit, but not so necessary
            end
            % wait until the name of the file in the buffer folder is
            % different from the previous read to make sure that the
            % observation frame is updated
            while(strcmp(char(strcat(obj.reader, '/', im_list(1))),obj.frameName)==1)
                pause(1);
                vid_ims = ls(obj.reader);
                im_list = cellstr(vid_ims(3:end,:));
            end
            % if the code is here, this means that the new file is ready
            % to be read
            fid = fopen(char(strcat(obj.reader, '/', im_list(1))), 'r'); % we try to read in case there is an operating system related lock on the file. This may happen if the camera is writing the
            % the file as we try to read it. if this application manages to open the file, than the file is locked by it and can be read.
            while (fid == -1) % if the camera is writing the file then apply the same trial based procedure until the write is over.
                if(read_try > 8)
                    disp('Input stream not available... Terminating...');
                    frame = [];
                    return;
                end
                disp(['Retrying from this one...   ' char(strcat(obj.reader, '/', im_list(1)))]);
                vid_ims = ls(obj.reader);
                im_list = cellstr(vid_ims(3:end,:));
                fid = fopen(char(strcat(obj.reader, '/', im_list(1))), 'r'); % try to lock the file once again
                read_try = read_try + 1;
                pause(0.5);
            end
            
            frame = imread(char(strcat(obj.reader, '/', im_list(1)))); % since the file is locked, it is safe to read it.
            fclose(fid); % since the file is read, the lock can be released
            
            obj.frameName = char(strcat(obj.reader, '/', im_list(1))); % record the name of the file for the check operation on the up-to-date information of the file
            obj.frameCount = obj.frameCount + 1; % increase the frame count
        else
            % if it is not a live operation, directly read the file in
            % order
            if(length(obj.im_list) >= obj.frameCount) % if there are still some images to be read
                obj.frameName = char(strcat(obj.reader, '/', obj.im_list(obj.frameCount)));
                frame = imread(obj.frameName);  % read the image directly, no need for the lock operations
                frameNameStr = strtok(char(obj.im_list(obj.frameCount)),'.'); % seperate the file type extension, get the date info
                frameDateVec = strread(frameNameStr,'%d','delimiter','_'); % get the date vector elements, in integer, corresponding to day, month etc.
                obj.frameDateVec = frameDateVec;
                
                if(mod(obj.frameCount-obj.frameRate,50*obj.frameRate)<=obj.frameRate)
                    frameDateVecStr = strread(frameNameStr,'%s','delimiter','_'); % get the date vector elements, in string, corresponding to day, month etc.
                    obj.frameRecordImageName = char(strcat(frameDateVecStr(1), '_',frameDateVecStr(2), '_',frameDateVecStr(3), ...
                                               '_',frameDateVecStr(4), '_',frameDateVecStr(5))); % the name of the image or the gif file with the directory
                end
                obj.frameCount = obj.frameCount + obj.frameRate;

            else
                frame = []; % if no more frames in the recording, directly terminate
                display('Done reading strored images... Terminating...');
            end
        end
        frame(~obj.imageMask) = 0;  % we eliminate the parts that are angularly behind the camera observation point
        %obj.predictionImage = frame;
        %obj.predictionImage(obj.imageMask) = 0;
        frame = imresize(frame,scale_factor);
    end

    function frame = readFrame_Simple()
        obj.frameCount = obj.frameCount + 1;
        frame = read(obj.reader,obj.frameCount);
    end

    function [contours, areas, centroids, bboxes, mask] = detectObjects(frame)
        
        % detect objects into a binary image
        mask = getBinaryImage(frame,maskMethod);
        mask = medfilt2(mask,[15 15]);
        
        % apply morphological operations to remove noise and fill in holes
        %mask = imopen(mask, strel('disk', 3));
        %mask = imclose(mask, strel('square', 30));
        %mask = imfill(mask, 'holes');
        
        contours = cv.findContours(mask); % contour extraction using openCV method
        
        %centroids = cellfun(@(x) mean(cell2mat(x(:)),1),contours,'UniformOutput',false);
        
        centroids = cellfun(@(x) (min(cell2mat(x(:)),[],1)+max(cell2mat(x(:)),[],1))/2,contours,'UniformOutput',false);
        
        bboxes = cellfun(@(x) [min(cell2mat(x(:)),[],1),max(cell2mat(x(:)),[],1)-min(cell2mat(x(:)),[],1)] ,contours,'UniformOutput',false);
        
        
        areas = cell2mat(cellfun(@(x) size(cell2mat(x(:)),1),contours,'UniformOutput',false))';
        
        % elminate the objects with contour smaller than 500.
        consider = areas > 300;
        
        centroids = centroids(consider);
        bboxes = bboxes(consider);
        areas = areas(consider);
        contours = contours(consider);
        
        centroids = vertcat(centroids{:});
        bboxes = vertcat(bboxes{:});
    end

    function [mask,cloud_mask] = getBinaryImage(frame, type, varargin)
        r_b = double(frame(:,:,1))./double(frame(:,:,3)); % R/B thresholding applied
        rgb = double(frame(:,:,1))+double(frame(:,:,2))+double(frame(:,:,3)); % to distinguish other objects like buildings from clouds, a test on brightness
        cloud_mask = r_b > 0.72 & rgb > 100;             %% brightness auto thing is needed because of this ..
        mask = double(frame(:,:,1));
        if(~isempty(varargin))
            sunP = varargin{1};
            g_b = double(frame(:,:,2))./double(frame(:,:,3));
            mask2 = g_b > 0.8;
            %mask2 = (r_b - g_b)>-0.01;
            cloud_mask(sunP) = mask2(sunP);
        end
        if(strcmp(type,'edge'))
            H = [1,1,1;1,-8,1;1,1,1];
            cloud_mask = imfilter(mask,H,'replicate');
        end
        cloud_mask = logical(cv.medianBlur(cloud_mask, 'KSize', 13));
        %mask = medfilt2(mask,[10 10]);
    end

    function [mask,cloud_mask] = getBinaryImage2(frame, type)
        
        r_b = double(frame(:,:,1))./double(frame(:,:,3)); % R/B thresholding applied
        cloud_mask = r_b > 0.72;
        mask = double(frame(:,:,1));
        
        cloudLevel = sum(sum(mask))/numel(mask);
        
        %         % if high cloud density, than we follow the sky segments..
        %         if(cloudLevel > 0.8)
        %             mask = r_b < 0.72;
        %             %mask = ones(size(mask)) - mask;
        %         end
        
        if(strcmp(type,'edge'))
            H = [1,1,1;1,-8,1;1,1,1];
            mask = imfilter(mask,H,'replicate');
            %mask = mask + mask2;
        end
        
    end

    function predictNewLocationsOfTracks()
        for i = 1:length(tracks)
            bbox = tracks(i).bbox;
            
            
            % predict the current location of the track
            [tracks(i).kalmanFilter,predictedCentroid,~] = tracks(i).kalmanFilter.Estimate();
            tracks(i).centroidFilter = filterAdd(predictedCentroid, tracks(i).centroidFilter);
            tracks(i).predCentroid = tracks(i).centroidFilter.value;
            
            switch trackOption
                case {'force', 'acceleration'}
                    tracks(i).predSpeed = [tracks(i).kalmanFilter.State(2),tracks(i).kalmanFilter.State(5)];
                case 'speed'
                    tracks(i).predSpeed = [tracks(i).kalmanFilter.State(2),tracks(i).kalmanFilter.State(4)];
                otherwise
            end
            
            
            predictedCentroid = tracks(i).centroidFilter.value;
            % shift the bounding box so that its center is at
            % the predicted location
            
            predictedCentroid = int32(predictedCentroid(1:2)) - (int32(bbox(3:4)) / 2);
            
            tracks(i).bbox = double([predictedCentroid(1:2), bbox(3:4)]);
        end
    end

    function [assignments, unassignedTracks, unassignedDetections] = ...
            detectionToTrackAssignment()
        nTracks = length(tracks);
        nDetections = size(centroids, 1);
        if(nTracks == 0)
            assignments = [];
            unassignedTracks = [];
            unassignedDetections = [1:nDetections];
            return;
        end
        
        % compute the cost of assigning each detection to each track
        cent_tracks = cat(1,tracks(:).centroid);
        match_mat1 = dot(cent_tracks,cent_tracks,2); % tracks
        match_mat2 = dot(centroids,centroids,2); % detections
        match_mat11 = repmat(match_mat1,1,length(match_mat2)); 
        match_mat22 = repmat(match_mat2',length(match_mat1),1);
        match_mat12 = cent_tracks*centroids';
        
        cost = match_mat11 - 2*match_mat12 + match_mat22;
        
        [matching,~] = Hungarian(cost);
        [assignments(:,1),assignments(:,2)] = ind2sub(size(matching),find(matching));
        unassignedTracks = setdiff([1:nTracks],assignments(:,1));
        unassignedDetections = setdiff([1:nDetections],assignments(:,2));
        
    end

    function updateBlindTracks()
        numAssignedTracks = size(assignments, 1);
        for i = 1:numAssignedTracks
            trackIdx = assignments(i, 1);
            % update track's age
            if(trackIdx > length(tracks))
                break;
            end
            tracks(trackIdx).age = tracks(trackIdx).age + 1;
            % update visibility
            tracks(trackIdx).consecutiveInvisibleCount = ...
                tracks(trackIdx).consecutiveInvisibleCount + 1;
            tracks(trackIdx).centroid = tracks(trackIdx).kalmanFilter.x';
        end
        for i = 1:length(unassignedTracks)
            ind = unassignedTracks(i);
            if(ind > length(tracks))
                break;
            end
            tracks(ind).age = tracks(ind).age + 1;
            tracks(ind).consecutiveInvisibleCount = ...
                tracks(ind).consecutiveInvisibleCount + 1;
            tracks(ind).centroid = tracks(trackIdx).kalmanFilter.x';
        end
    end

    function updateAssignedTracks()
        numAssignedTracks = size(assignments, 1);
        for i = 1:numAssignedTracks
            trackIdx = assignments(i, 1);
            detectionIdx = assignments(i, 2);
            centroid = centroids(detectionIdx, :);
            bbox = bboxes(detectionIdx, :);
            area = areas(detectionIdx, :);
            
            tracks(trackIdx).prevCentroid = tracks(trackIdx).centroid;
            %prevCentroid = tracks(trackIdx).kalmanFilter.State(1:2:3);
            % correct the estimate of the object's location
            % using the new detection
            
            %pos_t = (centroid - (frameSize(1:2)/2))*camParams.conversionFactor;
            
            pos_t = centroid;
            
            tracks(trackIdx).positionFilter = filterAdd(pos_t,tracks(trackIdx).positionFilter);
            pos_t = tracks(trackIdx).positionFilter.value;
            
            tracks(trackIdx).massFilter = filterAdd(area,tracks(trackIdx).massFilter);
            tracks(trackIdx).area = tracks(trackIdx).massFilter.value;
            
            %tracks(trackIdx).kalmanFilter = tracks(trackIdx).kalmanFilter.Update(pos_t',double(tracks(trackIdx).area));
            
            tracks(trackIdx).centroidFilter = filterAdd(centroid,tracks(trackIdx).centroidFilter);
            tracks(trackIdx).centroid = tracks(trackIdx).centroidFilter.value;
            
            tracks(trackIdx).obsSpeed = tracks(trackIdx).centroid - tracks(trackIdx).prevCentroid;
            
            tracks(trackIdx).kalmanFilter = tracks(trackIdx).kalmanFilter.Update(centroid',double(tracks(trackIdx).area));
            
            %             if(strcmp(trackOption,'areas'))
            %                 correct(tracks(trackIdx).kalmanFilter, [tracks(trackIdx).centroid,tracks(trackIdx).area]);
            %             else
            %                 correct(tracks(trackIdx).kalmanFilter, tracks(trackIdx).centroid);
            %             end
            
            % getting the speed estimate from the kalman...
            switch trackOption
                case {'force','acceleration'}
                    speed = [tracks(trackIdx).kalmanFilter.State(2), tracks(trackIdx).kalmanFilter.State(5)];
                case 'speed'
                    speed = [tracks(trackIdx).kalmanFilter.State(2), tracks(trackIdx).kalmanFilter.State(4)];
                otherwise
            end
            
            tracks(trackIdx).speedFilter = filterAdd(speed, tracks(trackIdx).speedFilter);
            speed = tracks(trackIdx).speedFilter.value;% * camParams.conversionFactor;
            
            tracks(trackIdx).xspeed = speed(2);
            tracks(trackIdx).yspeed = speed(1);
            
            % replace predicted bounding box with detected
            % bounding box
            tracks(trackIdx).bbox = bbox;
            
            % update track's age
            tracks(trackIdx).age = tracks(trackIdx).age + 1;
            
            % update visibility
            tracks(trackIdx).totalVisibleCount = ...
                tracks(trackIdx).totalVisibleCount + 1;
            tracks(trackIdx).consecutiveInvisibleCount = 0;
        end
    end

    function updateUnassignedTracks()
        for i = 1:length(unassignedTracks)
            ind = unassignedTracks(i);
            tracks(ind).age = tracks(ind).age + 1;
            tracks(ind).consecutiveInvisibleCount = ...
                tracks(ind).consecutiveInvisibleCount + 1;
        end
    end

    function deleteLostTracks()
        if isempty(tracks)
            return;
        end
        
        invisibleForTooLong = 900;
        ageThreshold = 1;
        
        % compute the fraction of the track's age for which it was visible
        ages = [tracks(:).age];
        totalVisibleCounts = [tracks(:).totalVisibleCount];
        visibility = totalVisibleCounts ./ ages;
        
        % find the indices of 'lost' tracks
        lostInds = (ages < ageThreshold & visibility < 0.8) | ...
            [tracks(:).consecutiveInvisibleCount] >= invisibleForTooLong;
        
        % delete lost tracks
        tracks = tracks(~lostInds);
    end

    function createNewTracks()
        centroids = centroids(unassignedDetections, :);
        bboxes = bboxes(unassignedDetections, :);
        areas = areas(unassignedDetections, :);
        
        for i = 1:size(centroids, 1)
            
            centroid = centroids(i,:);
            bbox = bboxes(i, :);
            
            area = areas(i);
            
            
            kalmanFilter = SimpleKF(centroid,1,[0.01,0.01],1e4,trackOption,double(area),1);
            % SimpleKF(init_x,sigma_p,sigma_q,sigma_r,type,m,time);
            
            frameSize = size(frame);
            %pos_t = (centroid - (frameSize(1:2)/2))*camParams.conversionFactor;
            
            % create a new track
            newTrack = struct(...
                'id', nextId, ...
                'bbox', bbox, ...
                'kalmanFilter', kalmanFilter, ...
                'age', 1, ...
                'totalVisibleCount', 1, ...
                'consecutiveInvisibleCount', 0, ...
                'xspeed', 0, ...
                'yspeed', 0, ...
                'predSpeed', [0,0], ...
                'obsSpeed', [0,0], ...
                'area', double(area), ...
                'speedFilter', filterInit(10,2), ...
                'positionFilter', filterInit(10,2), ...
                'centroidFilter', filterInit(1,2), ...
                'massFilter', filterInit(10,1), ...
                'centroid', centroid, ...
                'predCentroid', centroid, ...
                'prevCentroid', centroid);
            
            
            % add it to the array of tracks
            tracks(end + 1) = newTrack;
            
            % increment the next id
            nextId = nextId + 1;
        end
    end

    function displayTrackingResults()
        %         % convert the frame and the mask to uint8 RGB
        %         frame = im2uint8(frame);
        %
        %         mask = uint8(repmat(mask, [1, 1, 3])) .* 255;
        %
        minVisibleCount = 1;
        % draw on the mask
        frame = im2double(frame);
        frame = objectAnnotation(frame, 'CVArrow', arrows, 'r');
        %frame = objectAnnotation(frame, 'Rectangle', mesh_lut, 'y');
        if ~isempty(tracks)
            
            % noisy detections tend to result in short-lived tracks
            % only display tracks that have been visible for more than
            % a minimum number of frames.
            reliableTrackInds = ...
                [tracks(:).totalVisibleCount] > minVisibleCount;
            reliableTracks = tracks(reliableTrackInds);
            
            % display the objects. If an object has not been detected
            % in this frame, display its predicted bounding box.
            if ~isempty(reliableTracks)
                % get bounding boxes
                bboxes = cat(1, reliableTracks.bbox);
                
                % get ids
                %ids = int32([reliableTracks(:).id]);
                %xspeeds = [reliableTracks(:).xspeed];
                %yspeeds = [reliableTracks(:).yspeed];
                
                % create labels for objects indicating the ones for
                % which we display the predicted rather than the actual
                % location
                %labels = cellstr(int2str(ids'));
                %xlspeeds = cellstr(num2str(xspeeds',2));
                %ylspeeds = cellstr(num2str(yspeeds',2));
                
                %predictedTrackInds = ...
                %    [reliableTracks(:).consecutiveInvisibleCount] > 0;
                %isPredicted = cell(size(labels));
                %isPredicted(predictedTrackInds) = {' predicted'};
                %space = {', '};
                %labels = strcat(labels, ', x: ', xlspeeds, ', y: ', ylspeeds, isPredicted);
                predictedTrackInds = ...
                    [reliableTracks(:).consecutiveInvisibleCount] > 3;
                
                if (sum(predictedTrackInds)>0)
                    bboxColor = 'r';
                else
                    bboxColor = 'y';
                end
                % draw on the mask
                mask = objectAnnotation(mask, 'CVArrow', arrows, bboxColor);
                
                
                % draw on the frame
                frame = cv.drawContours(frame,contours,'Color',[255 255 0]);
                
            end
        end
        
        % draw on the mask
        
        %mask = objectAnnotation(mask, 'Rectangle', bboxes);
        
        % draw on the frame
        %frame = cv.drawContours(frame,contours,'Color',[255 255 0]);
        
        % display the mask and the frame
        obj.videoPlayer.Step(frame);
        obj.maskPlayer.Step(predictionImage);
        obj.predPlayer.Step(cloud_mask);
        %obj.videoPlayer.Step(frame);
    end

    function predictNewLocationsOfTracks_P()
        for i = 1:length(tracks2)
            bbox = tracks2(i).bbox;
            
            % predict the current location of the track
            predictedCentroid = predict(tracks2(i).kalmanFilter);
            
            tracks2(i).simpleKF = tracks2(i).simpleKF.Estimate();
            % shift the bounding box so that its center is at
            % the predicted location
            predictedCentroid = int32(predictedCentroid(1:2)) - bbox(3:4) / 2;
            
            tracks2(i).bbox = [predictedCentroid(1:2), bbox(3:4)];
        end
    end

    function [assignments, unassignedTracks, unassignedDetections] = ...
            detectionToTrackAssignment_P()
        
        nTracks = length(tracks2);
        nDetections = size(centroids2, 1);
        
        % compute the cost of assigning each detection to each track
        cost = zeros(nTracks, nDetections);
        for i = 1:nTracks
            if(strcmp(trackOption,'areas'))
                cost(i, :) = distance(tracks2(i).kalmanFilter, [centroids2,areas2]);
            else
                cost(i, :) = distance(tracks2(i).kalmanFilter, centroids2);
            end
        end
        
        % solve the assignment problem
        costOfNonAssignment = 1000000;
        [assignments, unassignedTracks, unassignedDetections] = ...
            assignDetectionsToTracks(cost, costOfNonAssignment);
    end

    function updateAssignedTracks_P()
        numAssignedTracks = size(assignments2, 1);
        frameSize = size(frame);
        for i = 1:numAssignedTracks
            trackIdx = assignments2(i, 1);
            detectionIdx = assignments2(i, 2);
            centroid = centroids2(detectionIdx, :);
            bbox = bboxes2(detectionIdx, :);
            area = areas2(detectionIdx, :);
            
            speed = centroid - tracks2(trackIdx).centroid;
            
            tracks2(trackIdx).xspeed = speed(2);
            tracks2(trackIdx).yspeed = speed(1);
            
            tracks2(trackIdx).centroid = centroid;
            
            tracks2(trackIdx).area = area;
            
            % replace predicted bounding box with detected
            % bounding box
            tracks2(trackIdx).bbox = bbox;
            
            % update track's age
            tracks2(trackIdx).age = tracks2(trackIdx).age + 1;
            
            % update visibility
            tracks2(trackIdx).totalVisibleCount = ...
                tracks2(trackIdx).totalVisibleCount + 1;
            tracks2(trackIdx).consecutiveInvisibleCount = 0;
        end
    end

    function updateUnassignedTracks_P()
        for i = 1:length(unassignedTracks2)
            ind = unassignedTracks2(i);
            tracks2(ind).age = tracks2(ind).age + 1;
            tracks2(ind).consecutiveInvisibleCount = ...
                tracks2(ind).consecutiveInvisibleCount + 1;
        end
    end

    function deleteLostTracks_P()
        if isempty(tracks2)
            return;
        end
        
        invisibleForTooLong = 900;
        ageThreshold = 8;
        
        % compute the fraction of the track's age for which it was visible
        ages = [tracks2(:).age];
        totalVisibleCounts = [tracks2(:).totalVisibleCount];
        visibility = totalVisibleCounts ./ ages;
        
        % find the indices of 'lost' tracks
        lostInds = (ages < ageThreshold & visibility < 0.8) | ...
            [tracks2(:).consecutiveInvisibleCount] >= invisibleForTooLong;
        
        % delete lost tracks
        tracks2 = tracks2(~lostInds);
    end

    function createNewTracks_P()
        centroids2 = centroids2(unassignedDetections2, :);
        bboxes2 = bboxes2(unassignedDetections2, :);
        areas2 = areas2(unassignedDetections2, :);
        
        for i = 1:size(centroids2, 1)
            
            centroid = centroids2(i,:);
            bbox = bboxes2(i, :);
            
            area = areas2(i);
            
            % create a Kalman filter object
            %kalmanFilter = configureKalmanFilter('ConstantVelocity', ...
            %    centroid, [10, 25], [100, 25], 100);
            if(strcmp(trackOption2,'areas'))
                state = [centroid(1), 0 , centroid(2), 0];
                kf_dims = 4;
                P = 1*eye(kf_dims);
                Q = 0.1 * eye(kf_dims);
                R = 20 * [1 0; 0 1]; % we have 2 measurements, so this should be in this shape.
                A = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]; %
                H = [1 0 0 0; 0 0 1 0];
                kalmanFilter = vision.KalmanFilter('StateTransitionModel', A, 'MeasurementModel', H, ...
                    'State', state, 'StateCovariance', P, 'ProcessNoise', Q, 'MeasurementNoise', R);
            elseif(strcmp(trackOption2,'pose_v'))
                kalmanFilter = configureKalmanFilter('ConstantVelocity', ...
                    centroid, [3, 3], [3, 3], 100);
            elseif(strcmp(trackOption2,'pose_acc'))
                kalmanFilter = configureKalmanFilter('ConstantAcceleration', ...
                    centroid, [10, 3, 5], [10, 1, 1], 100);
            end
            simpleKF = SimpleKF(10,1,20,double(area),1/29);
            
            frameSize = size(frame);
            pos_t = (centroid - (frameSize(1:2)/2))*camParams.conversionFactor;
            simpleKF = simpleKF.Update(pos_t(1),pos_t(2),double(area));
            
            % create a new track
            newTrack = struct(...
                'id', nextId, ...
                'bbox', bbox, ...
                'kalmanFilter', kalmanFilter, ...
                'simpleKF', simpleKF, ...
                'age', 1, ...
                'totalVisibleCount', 1, ...
                'consecutiveInvisibleCount', 0, ...
                'xspeed', 0, ...
                'yspeed', 0, ...
                'predSpeed', [0,0], ...
                'obsSpeed', [0,0], ...
                'area', double(area), ...
                'speedFilter', filterInit(10,2), ...
                'positionFilter', filterInit(10,2), ...
                'centroidFilter', filterInit(1,2), ...
                'massFilter', filterInit(10,1), ...
                'centroid', centroid, ...
                'predCentroid', centroid, ...
                'prevCentroid', centroid);
            %newTrack.rpixel = newTrack.xsensor/frameSize(2);
            %newTrack.rsensor = sqrt(newTrack.xsensor^2 + newTrack.ysensor^2);
            %newTrack.conversionFactor = sqrt((0.002/tan(45))^2 + (0.002^2 + 0.0013^2))/5000;
            % add it to the array of tracks
            tracks2(end + 1) = newTrack;
            
            % increment the next id
            nextId2 = nextId2 + 1;
        end
    end

    function frame = negativeFrame(frame)
        i_frame = ones(size(frame));
        frame = i_frame-frame;
    end
    
    function sunPosition = getSunPosition(frame,sunPattern)
        res = cv.matchTemplate(frame,sunPattern);
        if(min(min(res))>1e8)
            min(min(res))
        end
        [y,x]=ind2sub(size(res),find(res==min(min(res))));
        
        s_s = size(sunPattern);
        x_s = s_s(2);
        y_s = s_s(1);
        if(size(x)>(size(res,1)/3))
            sunPosition = [-1,-1];
            return;
        end
        sunPosition = [x,y]+([x_s,y_s]./2);
    end

    function frame = eliminateSun(frame,sunPosition)
        x = sunPosition(1);
        y = sunPosition(2);
        frame = cv.circle(frame,[x,y],80*scale_factor,'Color',[0,0,0],'Thickness',-1);
    end

    function sunPatch = getSunPatch(sunPosition)
        sunPatch = sqrt((obj.X-sunPosition(1)).^2+(obj.Y-sunPosition(2)).^2)<=200;
    end

    function pointReal = getImage2Real(m,calibModel,ecalibModel,cbh)
        % here we project a set of points given in m onto a plane that is cbh meters away from the camera the used convention in opencv is [column,row] kind of
        % represantation for pixel points, here however we change that to [row,column]
        pointTemp = cam2world([m(:,2),m(:,1)]',calibModel); % we get the corresponding direction on a unit-sphere at first.
        pointReal = zeros(size(pointTemp,2),2); % we crate the array that will contain the projected positions
        if(size(ecalibModel.R,1)==3)
            pointTemp = ecalibModel.R*pointTemp;
            pointReal = (pointTemp(1:2,:).*repmat((-cbh./pointTemp(3,:)),2,1))';
        else
            pointReal = ((ecalibModel.R*pointTemp(1:2,:)).*repmat((-cbh./pointTemp(3,:)),2,1))'; % here we map the point on the unit-sphere onto the plane
        end
        %for i = 1:size(pointTemp,2)
        %    % this part is a bit messy !!!! take care!!!/
        %    pointReal(i,:) = (ecalibModel.R*pointTemp(1:2,i))*(-cbh/pointTemp(3,i)); % here we map the point on the unit-sphere onto the plane
        %end
        %pointReal = mat2cell(pointReal,size(pointReal,1),size(pointReal,2));
    end

    function pointImage = getReal2Image(m,calibModel,ecalibModel,cbh)
        % here we project a set of points given in m back onto the image plane. The used convention in opencv is [column,row] kind of
        % represantation for pixel points, here however we change that to [row,column]
        if(size(ecalibModel.R,1)==3)
            pointTemp = world2cam(ecalibModel.Rinv*[m';-cbh*ones(1,size(m,1))],calibModel)';
        else
            pointTemp = world2cam([ecalibModel.R'*m';-cbh*ones(1,size(m,1))],calibModel)'; % we get the corresponding direction on a unit-sphere at first.
        end
        pointImage = [pointTemp(:,2),pointTemp(:,1)]-1;
    end

%% Generic filter implementation.
% c : length of the filter
% d : dimensions of the filter
% functions for initializing a filter, we can build a separate
% class out of this part...
    function filter = filterInit(c,d)
        filter = struct(...
            'capacity', c, ...
            'data', zeros(c,d), ...
            'cursor', 1, ...
            'numel', 0, ...
            'value', zeros(1,d));
    end

    function filter = filterAdd(newData, filter)
        
        filter.data(filter.cursor,:) = newData;
        filter.cursor = mod(filter.cursor + 1, filter.capacity)+1;
        if(filter.numel < filter.capacity)
            filter.numel = filter.numel + 1;
        end
        filter.value = sum(filter.data,1)/filter.numel;
    end
% end of filter implementation...

%% Camera filters and conversion factor implementation.
% xsensor   : size of the long edge of the camera sensor in m.
% ysensor   : size of the short edge of the camera sensor in m.
% rsensor   : size of the diagonal of the camera sensor in m.
% refsensor : size of the radius of the image on the sensor in m.
% rpixel    : size of a pixel in m.
% cbh       : estimated cloud base height for a frame in m.
% fov       : effective fov of the lens for the sensor in degrees.
% Camera parameter and cbh dependent conversion factor calculations. Can
% become a seperate file...
    function camParams = camParamsInit()
        
        xsensor = 0.0034;
        xefsensor = 0.0034;
        
        ysensor = 0.0019;
        yefsensor = 0.0019;
        
        rsensor = sqrt((xsensor/2)^2 + (ysensor/2)^2); %this depends on the camera sensor..
        refsensor = sqrt((xefsensor/2)^2 + (yefsensor/2)^2); %this depends on we utilize the camera sensor..
        
        frameRes = size(frame);
        
        rpixel = xsensor/frameRes(2);
        
        fov = 90 * 0.0174532925;
        
        cbh = 5000;
        
        timeConversionFactor = 1;
        lengthConversionFactor = cbh / sqrt(refsensor^2  + ((xefsensor/2)/tan(fov/2))^2);
        conversionFactor = rpixel * lengthConversionFactor / timeConversionFactor;
        
        camParams = struct(...
            'xsensor', xsensor, ...
            'ysensor', ysensor, ...
            'xefsensor', xefsensor, ...
            'yefsensor', yefsensor, ...
            'rsensor', rsensor, ...
            'refsensor', refsensor, ...
            'rpixel', rpixel, ...
            'lengthConversionFactor', lengthConversionFactor, ...
            'timeConversionFactor', timeConversionFactor, ...
            'conversionFactor', conversionFactor, ...
            'cbh', cbh, ...
            'xfov', fov);
    end
% end of camera parameter and cbh dependent part...

end
