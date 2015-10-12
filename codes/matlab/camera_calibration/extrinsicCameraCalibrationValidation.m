%% Cloud tracking using kalman filter.
%
% Multiple object tracking with kalman filter for could tracking:
%Inputs:
% image_folder: path to the input image files
% record: option indicating that a record should be taken or not
% verbose: option indication that the video should be seen or not
% calibModel: the intrinsic calibration model of the camera
% ecalibModel: the extrinsic calibration model of the camera
% imageMask: the mask indicating the usefull patch of the image
% sunPattern: the pattern for recognizing the sun location in the image
% live: whether we directly stream from the camera or from a folder
% track: whether we trackk or report using this function
% timeSeriesGHI: a database of the GHI data for some days
% timeSeriesPV: a database of the PV data for some days
% timeSeriesTemp: a database of temperature data for some days
%Outputs:
% recording_obs: Structure array that keeps the set of all observations.
% Kept mostly for evaluation and plotting purposes.
% recording_pred: Structure array that keeps the set of all predictions for
% all the timesteps. Kept mostly for evaluation and plotting purposes.
% db_objs: Recordingos of the observations.
% db_ass: Recordings of the observations to previously tracked objects.]

% B. Zeydan, 08. Oct. 2013

function [pix_err1, pix_err2, pix_errR1, pix_errR2, pix_errR3, pix_errR4, pix_errR5, pix_errR6 ]=extrinsicCameraCalibrationValidation(image_folder,calibModel,ecalibModel,model3D,imageMask,sunPattern,live,frameRate,R,cbh_given)
% create system objects used for 2reading video, detecting moving objects,
% and displaying the results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization of the key variables that are used as part of a model
Rinv = ecalibModel.Rinv;
sunPositionDetected = [];
sunPositionReal = [];
verbose = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the variables that determine the options and parameters applied on the
% tracking algorithm.

cbh = cbh_given; % the assumed cloud base height parameter that determines the projection performance of the algorithm.
sunPosition = [-1,-1]; % we initialize the position of sun, to be used throughout the code
deg2radf = pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% here we set-up the relevant objects.
obj = setupSystemObjects(image_folder,live,model3D,calibModel,ecalibModel,imageMask,cbh,frameRate);

if(isempty(obj.reader))
    disp('Camera does not respond. Terminating...');
    return;
end
disp_images = true;
iterator = 1;
% This is the main loop of the program
while obj.frameCount <= obj.finalFrame
    frame = readFrame();  % readFrame function returns the image, that is in order, from the buffer folder if this is not live, or the recording folder that is determined during the setup
    % This takes time, because it matches the sun pattern in the image, so
    % we run it not so frequently...
    sunPosition1 = getSunPosition(frame,sunPattern);
    %sunPositionDetected(:,iterator) = cam2world(sunPosition',obj.calibModel);
    sunPositionReal(:,iterator) = getSunPositionReal(obj.model3D);
    
    sunPosition = world2cam(Rinv*sunPositionReal(:,iterator),obj.calibModel);
    if (disp_images)
        figure;
        imshow(frame);
        hold on;
        plot(sunPosition(1),sunPosition(2),'Marker','o', 'MarkerSize',10, 'MarkerFaceColor','g');
        plot(sunPosition1(1),sunPosition1(2),'Marker','o', 'MarkerSize',10, 'MarkerFaceColor','b');
    end
    
%     sunPosition2 = getReal2Image(sunPositionReal(1:2,iterator)'/abs(sunPositionReal(3,iterator)),obj.calibModel,obj.ecalibModel,1);
%     pix_err1(iterator) = norm(sunPosition1-[sunPosition(2),sunPosition(1)]);
%     pix_err2(iterator) = norm(sunPosition1-sunPosition2);
%     
%     sunPositionCalculated = sunPositionReal(1:2,iterator)'*(-obj.cbh/sunPositionReal(3,iterator));
%     sunPositionDetected1 = R*(cam2world([sunPosition1(2),sunPosition1(1)]',obj.calibModel));
%     sunPositionDetected1 = (sunPositionDetected1(1:2).*repmat((-obj.cbh./sunPositionDetected1(3)),2,1))';
%     sunPositionDetected2 = getImage2Real(sunPosition1,obj.calibModel,obj.ecalibModel,obj.cbh);
%     
%     pix_errR1(iterator) = norm(sunPositionCalculated-sunPositionDetected1);
%     pix_errR2(iterator) = norm(sunPositionCalculated-sunPositionDetected2);
%     
%     sunPositionCalculated = sunPositionReal(1:2,iterator)'*(-(obj.cbh*1.01)/sunPositionReal(3,iterator));
%     
%     pix_errR3(iterator) = norm(sunPositionCalculated-sunPositionDetected1);
%     pix_errR4(iterator) = norm(sunPositionCalculated-sunPositionDetected2);
%     
%     sunPositionCalculated = sunPositionReal(1:2,iterator)'*(-(obj.cbh*0.99)/sunPositionReal(3,iterator));
%     
%     pix_errR5(iterator) = norm(sunPositionCalculated-sunPositionDetected1);
%     pix_errR6(iterator) = norm(sunPositionCalculated-sunPositionDetected2);
%     
%     frame = eliminateSun(frame,sunPosition);
%     frame = eliminateSunD(frame,sunPosition1);
%     frame = eliminateSunDD(frame,sunPosition2);
%     if(verbose)
%         displayTrackingResults();
%     end
    iterator = iterator + 1;
end

    function obj = setupSystemObjects(image_folder,live,model3D,calibModel,ecalibModel,imageMask,cbh,frameRate)
        % Initialize Video I/O
        % Create objects for reading a video from a file, drawing the tracked
        % objects in each frame, and playing the video.
        
        % create a video file reader
        try_ct = 0;
        obj.live = live;
        obj.model3D = model3D;
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
        % here some parameters to be used during the operation of the
        % algorithm are used.
        obj.gifMode = 0;    % if we want to export minutely recordings as a gif to the web page or not.
        obj.frameCount = 1; % the order of frames that are processed.
        obj.frameName = []; % the file name of the frame that is processed the last.
        obj.frameDateVec = [];
        obj.calibModel = calibModel; % the internal calibration model that is used by the algorithm
        obj.ecalibModel = ecalibModel; % the external calibration model that is used by the algorithm
        obj.cbh = cbh; % the assumed cloud base height that is used.
        obj.imageMask = imageMask;
        obj.frameRate = frameRate;
        if(verbose)
            obj.videoPlayer = VidPlayer([200, 500, 700, 400]); % for the real image
        end
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
            try
                frame = imresize(frame,[1800 1800]);
            catch    
                frame = cv.resize(frame, [1800 1800]);
            end
            frameNameStr = strtok(char(obj.im_list(obj.frameCount)),'.'); % seperate the file type extension, get the date info
            frameNameStr = strrep(frameNameStr, '_Debevec', '');
            frameDateVec = strread(frameNameStr,'%d','delimiter','_'); % get the date vector elements, in integer, corresponding to day, month etc.
            obj.frameDateVec = frameDateVec;
            obj.frameCount = obj.frameCount + obj.frameRate;
        else
            frame = []; % if no more frames in the recording, directly terminate
            display('Done reading strored images... Terminating...');
        end
        %frame(~obj.imageMask) = 0;  % we eliminate the parts that are angularly behind the camera observation point
    end
    
    function [contours, contoursProj,areas, centroids, centroidsProj, bboxes, bboxesProj, mask] = detectObjects(frame)
        % detect objects into a binary image
        mask = getBinaryImage(frame,maskMethod);
        mask = medfilt2(mask,[5 5]);
        % apply morphological operations to remove noise and fill in holes
        %mask = imopen(mask, strel('disk', 3));
        %mask = imclose(mask, strel('square', 30));
        %mask = imfill(mask, 'holes');
        contours = cv.findContours(mask,'Mode','List'); % contour extraction using openCV method
        % areas of objects are approximately computed here using several methods
        %areas = cell2mat(cellfun(@(x) size(cell2mat(x(:)),1),contours,'UniformOutput',false))';
        %areas = cell2mat(cellfun(@(x) prod(max(cell2mat(x(:)),[],1)-min(cell2mat(x(:)),[],1)) ,contours,'UniformOutput',false))';
        areas = cell2mat(cellfun(@(x) cv.contourArea(x(:)) ,contours,'UniformOutput',false))';
        % elminate the objects with contour smaller than 350.
        consider = (areas > 250) & (areas < 100500); % if the area of an object is less than 350 pixels, we discard that
        
        % take only the objects that we want to consider
        contours = contours(consider);
        areas = areas(consider);
        
        contoursProj = cellfun(@(x) num2cell(getImage2Real(cell2mat(x(:)),obj.calibModel,obj.ecalibModel,obj.cbh),2)',contours,'UniformOutput',false);
        
        % different methods for determining the point that represent the
        % cloud objets first method uses the mean of the contours the
        % second one uses the mid point of the bounding rectangle
        %centroids = cellfun(@(x) mean(cell2mat(x(:)),1),contours,'UniformOutput',false);
        centroids = cellfun(@(x) (min(cell2mat(x(:)),[],1)+max(cell2mat(x(:)),[],1))/2,contours,'UniformOutput',false);
        centroidsProj = cellfun(@(x) (min(cell2mat(x(:)),[],1)+max(cell2mat(x(:)),[],1))/2,contoursProj,'UniformOutput',false);
        
        bboxes = cellfun(@(x) [min(cell2mat(x(:)),[],1),max(cell2mat(x(:)),[],1)-min(cell2mat(x(:)),[],1)] ,contours,'UniformOutput',false);
        bboxesProj = cellfun(@(x) [min(cell2mat(x(:)),[],1),max(cell2mat(x(:)),[],1)-min(cell2mat(x(:)),[],1)] ,contoursProj,'UniformOutput',false);
        
        centroids = vertcat(centroids{:});
        centroidsProj = vertcat(centroidsProj{:});
        bboxes = vertcat(bboxes{:}); % put all the boxes in a 2D array, so that the annotation algorithm can display them
        bboxesProj = vertcat(bboxesProj{:});
    end

    function sunPosition = getSunPosition(frame,sunPattern)        
        res = cv.matchTemplate(frame,sunPattern);
        matchQ = min(min(res));
        [y,x]=ind2sub(size(res),find(res==matchQ));
        s_s = size(sunPattern);
        y_s = s_s(1);
        x_s = s_s(2);
        %if(size(x)>(size(res,1)/3))
        %    sunPosition = [-1,-1];
        %    return;
        %end
        if(length(x)>1)
            sunPosition = [-1,-1];
            return;
        end
        sunPosition = [x,y]+([x_s,y_s]./2);
        
    end

    function frame = eliminateSun(frame,sunPosition)
        x = sunPosition(1);
        y = sunPosition(2);
        if(x<0 || y<0)
            return;
        end
        frame = cv.circle(frame,[y,x],20,'Color',[255,0,0],'Thickness',-1);
    end
    function frame = eliminateSunD(frame,sunPosition)
        x = sunPosition(1);
        y = sunPosition(2);
        if(x<0 || y<0)
            return;
        end
        frame = cv.circle(frame,[x,y],20,'Color',[0,0,255],'Thickness',-1);
    end
    function frame = eliminateSunDD(frame,sunPosition)
        x = sunPosition(1);
        y = sunPosition(2);
        if(x<0 || y<0)
            return;
        end
        frame = cv.circle(frame,[x,y],20,'Color',[0,255,255],'Thickness',-1);
    end

    function sunPositionReal = getSunPositionReal(model3D)
        DN = datenum(obj.frameDateVec(1:6)');
        Time = pvl_maketimestruct(DN, model3D.UTC);
        %[SunAz, ~, ApparentSunEl, ~]=pvl_ephemeris(Time, model3D.Location);
        [SunAz, ~, ApparentSunEl]=pvl_spa(Time, model3D.Location);
        zenith = 90-ApparentSunEl;
        azimuth = SunAz;
        [x,y,z] = sph2cart((90-azimuth)*deg2radf,(90-zenith)*deg2radf,1);
        sunPositionReal = [x;y;-z];
    end

    function pointReal = getImage2Real(m,calibModel,ecalibModel,cbh)
        % here we project a set of points given in m onto a plane that is cbh meters away from the camera the used convention in opencv is [column,row] kind of
        % represantation for pixel points, here however we change that to [row,column]
        pointTemp = cam2world([m(:,2),m(:,1)]',calibModel); % we get the corresponding direction on a unit-sphere at first.
        pointReal = zeros(size(pointTemp,2),2); % we crate the array that will contain the projected positions
        for i = 1:size(pointTemp,2)
            % this part is a bit messy !!!! take care!!!/
            pointReal(i,:) = (ecalibModel.R*pointTemp(1:2,i))*(-cbh/pointTemp(3,i)); % here we map the point on the unit-sphere onto the plane
        end
        %pointReal = mat2cell(pointReal,size(pointReal,1),size(pointReal,2));
    end

    function pointImage = getReal2Image(m,calibModel,ecalibModel,cbh)
        % here we project a set of points given in m back onto the image plane. The used convention in opencv is [column,row] kind of
        % represantation for pixel points, here however we change that to [row,column]
        pointTemp = world2cam([ecalibModel.R'*m';-cbh*ones(1,size(m,1))],calibModel)'; % we get the corresponding direction on a unit-sphere at first.
        pointImage = [pointTemp(:,2),pointTemp(:,1)];
    end
    
    function displayTrackingResults()
        obj.videoPlayer.Step(frame);
    end
end