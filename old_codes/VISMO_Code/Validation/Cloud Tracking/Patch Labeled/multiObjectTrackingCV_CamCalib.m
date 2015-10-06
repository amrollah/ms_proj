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

function [recording_obs,recording_pred,db_objs,db_ass] = multiObjectTrackingCV_CamCalib(image_folder,record,verbose,calibModel,ecalibModel,model3D,imageMask,sunPattern,modelImage,live,track,timeSeriesGHI,timeSeriesPV,timeSeriesTemp)
% create system objects used for 2reading video, detecting moving objects,
% and displaying the results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization of some input variables in case they are not provided to
% the method..
if (nargin < 2)
    record = 0;
end

if(nargin < 3)
    verbose = 1;
end

if(nargin < 4)
    finalFrame = 10000;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization of the key variables that are used as part of a model

tracks = initializeTracks(); % create an empty array of tracks
nextId = 1; % ID of the next track
centroids = []; % the array that keeps the centroid positions of tracks on the image plane
centroidsProj = []; % the array that keeps the centroid positions of tracks on the projected (real) plane
bboxes = []; % the array that keeps the bounding boxes of the detected cloud trails in the image plane
bboxesProj = []; % the array that keeps the bounding boxes of the detected cloud trails in the projected (real) plane
areas = []; % the array that keeps the areas of bounding boxes
mask = []; % the element that keeps the binary image that contains clouds only
predCloudContours = [];
predShadowContours = [];
% these are the variables that are used as a buffer for the validation of
% the kalman filter for tracking
if(record && (predictionCycle > updateCycle))
    tracks2 = initializeTracks(); % create an empty array of tracks
    nextId2 = 1; % ID of the next track
    trackOption2 = 'pose_v'; % the option that is applied for the kalman filter, for other possibilities please look inside the object implementation.
    centroids2 = []; % the array that keeps the centroid positions of tracks on the image plane
    centroidsProj2 = []; % the array that keeps the centroid positions of tracks on the projected (real) plane
    bboxes2 = []; % the array that keeps the bounding boxes of the detected cloud trails in the image plane
    bboxesProj2 = []; % the array that keeps the bounding boxes of the detected cloud trails in the projected (real) plane
    areas2 = []; % the array that keeps the areas of bounding boxes
    mask2 = []; % the element that keeps the binary image that contains clouds only
end

% in case we go with the validation mode and we are recording the
% performance of the algorithm, the pstate variables of all the steps are
% recorded in these variables
if(record)
    recording_pred = cell(finalFrame,1);
    recording_obs = cell(finalFrame,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the variables that determine the options and parameters applied on the
% tracking algorithm.

cbh = 300; % the assumed cloud base height parameter that determines the projection performance of the algorithm.
ModPlayerOption = 'map'; % the representation of clouds on the rea-sized model is done using this opttion. rectangle is another possibility.
maskMethod = 'body'; % identification type of clouds is done by this variable. either they are identified by their body or the contours (edges)
trackOption = 'speed'; % the option that is applied for the kalman filter, for other possibilities please look inside the object implementation.
predictionCycle = 8; % the algorithm has the sense of periods for the following reason: we want to validate that the tracking algorithm works
% properly for this reason, we update the kalman filters for several steps, and then we try to test the predictions of
% kalman filter not updating with the observations but instead going with the predictions. This is the period length.
% For the first 'updateCycle' of this period, the filter are updated and for the rest they are not
updateCycle = 8; % the number of cycles where the kalman filter is updated by observations
sunEliminationCycle = 100; % every this many frames, we calculate the sun's position once again. 
sunPosition = [-1,-1]; % we initialize the position of sun, to be used throughout the code
deg2radf = pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% here we set-up the relevant objects.
obj = setupSystemObjects(image_folder,live,verbose,calibModel,ecalibModel,model3D,imageMask,modelImage,cbh,ModPlayerOption,track,timeSeriesGHI,timeSeriesPV,timeSeriesTemp);

if(isempty(obj.reader))
    disp('Camera does not respond. Terminating...');
    return;
end

db_objs = [];
db_ass = [];


% This is the main loop of the program
while obj.frameCount < obj.finalFrame
    
    frame = readFrame();  % readFrame function returns the image, that is in order, from the buffer folder if this is not live, or the recording folder that is determined during the setup
   
    % This takes time, because it matches the sun pattern in the image, so
    % we run it not so frequently...
    if((sum(sunPosition)<0) || mod(obj.frameCount,sunEliminationCycle)>=(sunEliminationCycle-obj.frameRate))
        sunPosition = getSunPosition(frame,sunPattern);
        getShadowDisplacement(model3D)
    end
    
    if(isempty(frame))    % if the folder does not exist or the frame is not correct, we immediately terminate
        disp('Camera does not respond. Terminating...');
        close all;
        return;
    end
    
    camParams = camParamsInit(); % this is a first implementation of the camera parameters that is used for camera frame to real world frame. not used any more
    if(track) % a special option for this file, if the user wants to track, that is the option, else if he wants to dump the collected data into a web server
        
        if((mod(obj.frameCount,predictionCycle) < updateCycle) || (obj.frameCount < 3)) % this is the part that controls the kalman update-predict cycles. If this
                                                                                        % statement is correct than the set of kalman filters are updated.
            frame = eliminateSun(frame,sunPosition);
            [contours, contoursProj, areas, centroids, centroidsProj, bboxes, bboxesProj, mask] = detectObjects(frame); % this function takes in the RGB frame and :
                                                                                                                        % - runs the cloud decision algorithm to get a binary image
                                                                                                                        % - runs a clustering algorithm followed by a contour extraction
                                                                                                                        % - runs simple heuristic methods on the contour pixels to extract
                                                                                                                        %   - the representative points for the model based tracking
                                                                                                                        %   - the projected representative point for physical interpretation
                                                                                                                        %   - bounding boxes for the objects found
                                                                                                                        % - returns all the results along with the binary image to be used for
                                                                                                                        % displaying
            
            if(record)
                db_objs{obj.frameCount} = struct('contours',contours,...  % Here we save the current state of the algorithm in order to be able to run the monitor the algorithm better
                    'areas',areas,...                                     % This corresponds to the "obeservations" part of the system from a Kalman perspective
                    'centroids',centroids,...
                    'bboxes',bboxes);
            end
            
            predictNewLocationsOfTracks();    % in this function we predict the one-step future of the existing model calling the predict method of
                                              % existing Kalman filters.
                                              
            predTracks = predictNewLocationsOfTracksIntoFuture(20);
            predFrame = getPredictedImage(predTracks);
            [predCloudContours,predShadowContours] = getPredictedScenario(predTracks);
            
            
            if(record && (predictionCycle >= updateCycle))
                recording_pred{obj.frameCount} = tracks;  % this is the database where we record the predictions of the model before the updates.
            end
            
            [assignments, unassignedTracks, unassignedDetections] = ...
                detectionToTrackAssignment(); % for every frame that is processed, there is a set of detections. Also for the frames previous to the current one, there are the tracks.
            % This function establishes the association between the observations and tracks.
            
            if(record)
                db_ass{obj.frameCount} = struct('assignments', assignments,... % this is the database file that keeps track of the assignments.
                    'unassignedTracks', unassignedTracks,...
                    'unassignedDetections', unassignedDetections);
            end
            
            updateAssignedTracks(); % this is the function that updates the tracks according to assignments, it also updates some statistical variables related to update rate
            updateUnassignedTracks(); % this is the funstion that updates tracks corresponding to unassigned tracks, it just updates statistical variables related to update rate
            deleteLostTracks(); % here we delete tracks that remained unassigned for some number of steps
            createNewTracks(); % here we create new tracks corresponding to unassigned observations
        else
            % this part of the code runs only for the cycles that do not
            % update kalman filters
            mask = getBinaryImage(frame,maskMethod); % here, for displaying purposes we get the binary image of the current frame
            mask = medfilt2(mask,[10 10]); % here we apply the very same filtering on the binary
            
            predictNewLocationsOfTracks(); % here we call prediction on kalman filters
            updateBlindTracks(); % here we update statistical paramaters of all the tracks, since they are not updated
            %deleteLostTracks();
        end
        
        % for recording purposes...
        if(record)
            if(obj.frameCount == 1)
                recording_pred{obj.frameCount} = tracks;
            end
            recording_obs{obj.frameCount} = tracks; % here we record tracks after update
        end
    end
    if(verbose)
        displayTrackingResults(); % here we display the status and results
    end
end

    function obj = setupSystemObjects(image_folder,live,verbose,calibModel,ecalibModel,model3D,imageMask,modelImage,cbh,ModPlayerOption,track,timeSeriesGHI,timeSeriesPV,timeSeriesTemp)
        % Initialize Video I/O
        % Create objects for reading a video from a file, drawing the tracked
        % objects in each frame, and playing the video.
        
        % create a video file reader
        try_ct = 0;
        obj.live = live;
        
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
        obj.calibModel = calibModel; % the internal calibration model that is used by the algorithm
        obj.ecalibModel = ecalibModel; % the external calibration model that is used by the algorithm
        obj.cbh = cbh; % the assumed cloud base height that is used.
        obj.ModPlayerOption = ModPlayerOption; % the way clouds represented in the model player (can be cycle or rectangle)
        obj.track = track; % whether we are in the tracking mode or display mode. can be 0 or 1
        obj.timeSeriesGHI = timeSeriesGHI; % the timeseries object corresponding to the GHI measurements
        obj.timeSeriesPV = timeSeriesPV; % the timeseries object corresponding to the PV measurements
        obj.timeSeriesTemp = timeSeriesTemp; % the timeseries object corresponding to the temperature measurements
        obj.irradianceValue = 0; % the variable that keeps the irradiance value corresponding to the current frame
        obj.imageMask = imageMask;
        obj.frameRate = 4;
        if(obj.track) % if the tracking mode is active initialize the needed display objects
            % create three video players, one to display the video,
            % one to display the foreground mask and objects, and one
            % to display the model.
            if(verbose)
                obj.videoPlayer = VidPlayer([20, 500, 700, 400]); % for the real image
                obj.maskPlayer = VidPlayer([740, 500, 700, 400]); % for the masked cloud-decision image
                %obj.modPlayer = ModPlayer([740, 100, 700, 400],obj.ModPlayerOption); % for the model in real dimensions
                obj.modPlayer = [];
                obj.predPlayer = VidPlayer([740, 100, 700, 400]); % for the model in real dimensions
                obj.mapPlayer = ModPlayer([20, 100, 700, 400],obj.ModPlayerOption,modelImage,model3D.all); % for the model in real dimensions
            end
        else % othervise initialize the other objects
            % create two video players, one to display the video,
            % and one to display the irradiance values
            if(verbose)
                obj.videoPlayer = VidPlayer([20, 500, 700, 300]); % for the real image
                obj.plotPlayerGHI = ModPlayer([740, 500, 700, 300],'plot'); % for the GHI-value plot
                obj.plotPlayerPV = ModPlayer([20, 50, 700, 300],'plot'); % for the PV-value plot
                obj.plotPlayerTemp = ModPlayer([740, 50, 700, 300],'plot'); % for the temperature-value plot
            end
        end
    end

    function tracks = initializeTracks()
        % create an empty array of tracks
        tracks = struct(...
            'id', {}, ...   % id of the track
            'bbox', {}, ... % the bounding box corresponding to object, on the image axis
            'bboxProj', {}, ... % the bounding box corresponding to object, projected onto real-world axis
            'kalmanFilter', {}, ... % kalman filter corresponding to this object
            'kalmanFilterProj', {}, ... % kalman filter corresponding to this object, projected onto real-world axis
            'age', {}, ... % how old the track is
            'totalVisibleCount', {}, ... % for how many cycles do we get an observation on this track
            'consecutiveInvisibleCount', {},... % if we cannot observe it, for how many previous steps we could not
            'xspeed', {}, ... % a variable for keeping one directional speed in image plane
            'yspeed', {}, ... % a variable for keeping one directional speed in image plane
            'predSpeed', {}, ... % a variable keeping the 2D predicted speed in image plane
            'obsSpeed', {}, ... % a variable keeping the 2D observed speed in image plane
            'area', {}, ... % a variable keeping the area of the tracked object
            'contour', {}, ... % a variable keeping the relevant contour points on the image plane
            'contourProj', {}, ... % a variable keeping the relevant contour points on the real world coordinates
            'speedFilter', {}, ... % a filter for the observed speed on the image plane
            'positionFilter', {}, ... % a filter for the position on the image plane
            'centroidFilter', {}, ... % a filter for the centroid on the image plane
            'centroidProjFilter', {}, ... % a filter for the centroid on the real plane
            'massFilter', {}, ... % a filter on the predicted mass (not used)
            'centroid', {}, ... % the representative position of the center of the object on the image plane
            'centroidProj', {}, ... % the representative position of the center of the object on the projected plane
            'predCentroid', {}, ... % prediction on the representative position of the center of the object on the image plane
            'predCentroidProj', {}, ... % prediction on the representative position of the center of the object on the projected plane
            'prevCentroid', {}, ... % previous centroid position on the image plane
            'prevCentroidProj', {}); % previous centroid position on the projected plane
        
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
        
        if(~obj.track)
            % prepare the elements for the required date strings
            obj.frameName = char(strcat(obj.reader, '/', im_list(obj.frameCount))); % get the file name for the last frame
            frameNameStr = strtok(char(im_list(obj.frameCount)),'.'); % seperate the file type extension, get the date info
            frameDateVec = strread(frameNameStr,'%d','delimiter','_'); % get the date vector elements, in integer, corresponding to day, month etc.
            frameDateVecStr= strread(frameNameStr,'%s','delimiter','_'); % get the date vector elements, in string, corresponding to day, month etc.
            frameDateVec = frameDateVec(1:6); % eliminate the unnecessary part of the frame vector
            
            % if does not exist, build the correct directory in the web server
            if(exist(char(strcat(frameDateVecStr(1), '/', frameDateVecStr(2), '/', frameDateVecStr(3),'/images/')))~=7)
                mkdir(char(strcat(frameDateVecStr(1), '/', frameDateVecStr(2), '/', frameDateVecStr(3),'/images/'))); % make the directory with all the sub directories
            end
            
            % build required date strings
            if(mod(obj.frameCount,20)==2)
                obj.frameRecordImageName = char(strcat(frameDateVecStr(1), '/', frameDateVecStr(2), '/', frameDateVecStr(3), '/images/',...
                    frameDateVecStr(1), '_',frameDateVecStr(2), '_',frameDateVecStr(3), '_',frameDateVecStr(4), '_',...
                    frameDateVecStr(5))); % the name of the image or the gif file with the directory
                obj.frameRecordImagePageName = char(strcat('images/',frameDateVecStr(1), '_',frameDateVecStr(2),'_',...
                    frameDateVecStr(3), '_',frameDateVecStr(4), '_',frameDateVecStr(5))); % the name of the image or the gif file without the directory
                obj.frameRecordPageName = char(strcat(frameDateVecStr(1), '/', frameDateVecStr(2), '/', frameDateVecStr(3), '/',...
                    frameDateVecStr(1), '_',frameDateVecStr(2), '_',frameDateVecStr(3), '_',frameDateVecStr(4), '_',...
                    frameDateVecStr(5))); % the name of the  web page file with the directory
            end
            
            % get the frame date in string format to query the needed period from timeseries objects
            frameDateStr = datestr(frameDateVec');
            obj.frameDateStr = frameDateStr;
            minuteLag = 20; % go back for this much of minutes
            frameDateStr2 = datestr(datenum(frameDateStr)-(minuteLag/(60*24)));
            
            % read the irradiance value between the current frame date and minuteLag mins before
            irradianceDB = obj.timeSeriesGHI.getsampleusingtime(frameDateStr2,frameDateStr); % query the values
            obj.irradianceValues = irradianceDB.Data*100; % scale the values
            obj.irradianceTimes = ( (irradianceDB.Time + datenum(irradianceDB.TimeInfo.StartDate)) - datenum(frameDateStr) )*60*24; % get the time in terms of difference between the current time
            if(isempty(obj.irradianceValues)) % if not found the associated entries, set to zero
                obj.irradianceValues = 0;
                obj.irradianceTimes = 0;
            end
            
            % read the power value between the current frame date and minuteLag mins before.
            powerDB = obj.timeSeriesPV.getsampleusingtime(frameDateStr2,frameDateStr); % query the values
            obj.powerValues = powerDB.Data; % get the values into buffer
            obj.powerTimes = ( (powerDB.Time + datenum(powerDB.TimeInfo.StartDate)) - datenum(frameDateStr) )*60*24; % get the time in terms of difference between the current time
            if(isempty(obj.powerValues)) % if not found the associated entries, set to zero
                obj.powerValues = 0;
                obj.powerTimes = 0;
            end
            
            % read the temperature value between the current frame date and minuteLag mins before.
            tempDB = obj.timeSeriesTemp.getsampleusingtime(frameDateStr2,frameDateStr); % query the values
            obj.tempValues = tempDB.Data; % get the values into buffer
            obj.tempTimes = ( (tempDB.Time + datenum(tempDB.TimeInfo.StartDate)) - datenum(frameDateStr) )*60*24; % get the time in terms of difference between the current time
            if(isempty(obj.tempValues)) % if not found the associated entries, set to zero
                obj.tempValues = 0;
                obj.tempTimes = 0;
            end
        end
        frame(~obj.imageMask) = 0;  % we eliminate the parts that are angularly behind the camera observation point
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

    function mask = getBinaryImage(frame, type)
        r_b = double(frame(:,:,1))./double(frame(:,:,3)); % R/B thresholding applied
        rgb = double(frame(:,:,1))+double(frame(:,:,2))+double(frame(:,:,3)); % to distinguish other objects like buildings from clouds, a test on brightness
        mask = r_b > 0.72 & rgb > 100;             %% brightness auto thing is needed because of this ..
        if(strcmp(type,'edge'))
            H = [1,1,1;1,-8,1;1,1,1];
            mask = imfilter(mask,H,'replicate');
        end        
    end

    function sunPosition = getSunPosition(frame,sunPattern)
        res = cv.matchTemplate(frame,sunPattern);
        [x,y]=ind2sub(size(res),find(res==min(min(res))));
        sunPosition = [x,y];
    end

    function frame = eliminateSun(frame,sunPosition)
        x = sunPosition(1);
        y = sunPosition(2);
        frame = cv.circle(frame,[y+3,x+1],130,'Color',[0,0,0],'Thickness',-1);
    end
    
    function getShadowDisplacement(model3D)
        DN = datenum(obj.frameDateVec(1:6)');
        Time = pvl_maketimestruct(DN, model3D.UTC);
        [SunAz, ~, ApparentSunEl, ~]=pvl_ephemeris(Time, model3D.Location);
        zenith = 90-ApparentSunEl;
        azimuth = SunAz;
        [x,y,~] = sph2cart((90-azimuth)*deg2radf,(90-zenith)*deg2radf,cbh/cos(zenith*deg2radf));
        obj.shadowDisp = [-x,-y];
    end

    function tracksPred = predictNewLocationsOfTracksIntoFuture(timeSteps)
        tracksPred = tracks;
        for predNum = 1:timeSteps
            for i = 1:length(tracksPred)
                % to be able to used them for image annotation and sanity
                % checking, we provide the bounding boxes in an array format
                % where one bounding box is on one row of the matrix
                bbox = tracksPred(i).bbox;
                bboxProj = tracksPred(i).bboxProj;
                
                % predict the current location of the track
                [tracksPred(i).kalmanFilter,predictedCentroid,~] = tracksPred(i).kalmanFilter.Estimate();
                tracksPred(i).centroidFilter = filterAdd(predictedCentroid, tracksPred(i).centroidFilter);  % the extra filters are not operating right now, they are put in case needed.
                tracksPred(i).predCentroid = tracksPred(i).centroidFilter.value;
                
                % projected values from the kalman filter.
                [tracksPred(i).kalmanFilterProj,predictedCentroidProj,~] = tracksPred(i).kalmanFilterProj.Estimate();
                tracksPred(i).centroidProjFilter = filterAdd(predictedCentroidProj, tracksPred(i).centroidProjFilter);  % the extra filters are not operating right now, they are put in case needed.
                tracksPred(i).predCentroidProj = tracksPred(i).centroidProjFilter.value;
                
                % here we update the state variables of the current track
                % according to different types of the kalman filter
                switch trackOption
                    case {'force', 'acceleration'}
                        tracksPred(i).predSpeed = [tracksPred(i).kalmanFilter.State(2),tracksPred(i).kalmanFilter.State(5)];
                    case 'speed'
                        tracksPred(i).predSpeed = [tracksPred(i).kalmanFilter.State(2),tracksPred(i).kalmanFilter.State(4)];
                    otherwise
                end
                
                % here to update bounding boxes, we get a local copy on the
                % location of tracks
                predictedCentroid = tracksPred(i).predCentroid;
                
                % shift the bounding box so that its center is at
                % the predicted location
                predictedCentroid = int32(predictedCentroid(1:2)) - (int32(bbox(3:4)) / 2);
                
                % reform  the bounding boxes with their new locations
                tracksPred(i).bbox = double([predictedCentroid(1:2), bbox(3:4)]);
                
                % here to update bounding boxes (on the projected plane), we get a local copy on the
                % location of tracks
                predictedCentroidProj = tracksPred(i).predCentroidProj;
                
                % shift the bounding box so that its center is at
                % the predicted location
                switch obj.ModPlayerOption
                    case 'rectangle'
                        predictedCentroidProj = int32(predictedCentroidProj(1:2)) - (int32(bboxProj(3:4)) / 2);
                        predictedCentroidProj(2) = predictedCentroidProj(2)-bboxProj(4);
                        tracksPred(i).bboxProj = double([predictedCentroidProj(1:2), bboxProj(3:4)]);
                    case 'circle'
                        tracksPred(i).bboxProj = double([predictedCentroidProj(1:2), bboxProj(3)]);
                    otherwise
                end
            end
        end
    end

    function predFrame = getPredictedImage(tracksPred)
        predFrame = zeros(size(frame));
        predContours=[];
        for i = 1:length(tracksPred)
            predContours{1,i} = num2cell(getReal2Image(tracksPred(i).contourProj+repmat(tracksPred(i).predCentroidProj,size(tracksPred(i).contourProj,1),1),obj.calibModel,obj.ecalibModel,cbh),2)'; 
        end
        if(length(tracksPred)~=0)
            predFrame = cv.drawContours(predFrame,predContours,'Thickness',-1,'LineType',4);
        end
    end

    function [predCloudContours,predShadowContours] = getPredictedScenario(tracksPred)
        predCloudContours=[];
        predShadowContours=[];
        for i = 1:length(tracksPred)
            predCloudContours{1,i} = num2cell([tracksPred(i).contourProj+repmat(tracksPred(i).predCentroidProj,size(tracksPred(i).contourProj,1),1),cbh*ones(size(tracksPred(i).contourProj,1),1)],2)'; 
            predShadowContours{1,i} = num2cell([tracksPred(i).contourProj+repmat(tracksPred(i).predCentroidProj+obj.shadowDisp,size(tracksPred(i).contourProj,1),1),zeros(size(tracksPred(i).contourProj,1),1)],2)';
        end
    end

    function predictNewLocationsOfTracks()
        for i = 1:length(tracks)
            % to be able to used them for image annotation and sanity
            % checking, we provide the bounding boxes in an array format
            % where one bounding box is on one row of the matrix
            bbox = tracks(i).bbox;
            bboxProj = tracks(i).bboxProj;
            
            % predict the current location of the track
            [tracks(i).kalmanFilter,predictedCentroid,~] = tracks(i).kalmanFilter.Estimate();
            tracks(i).centroidFilter = filterAdd(predictedCentroid, tracks(i).centroidFilter);  % the extra filters are not operating right now, they are put in case needed.
            tracks(i).predCentroid = tracks(i).centroidFilter.value;
            
            % projected values from the kalman filter.
            [tracks(i).kalmanFilterProj,predictedCentroidProj,~] = tracks(i).kalmanFilterProj.Estimate();
            tracks(i).centroidProjFilter = filterAdd(predictedCentroidProj, tracks(i).centroidProjFilter);  % the extra filters are not operating right now, they are put in case needed.
            tracks(i).predCentroidProj = tracks(i).centroidProjFilter.value;
            
            %tracks(i).predCentroid = getReal2Image(tracks(i).predCentroidProj,obj.calibModel,obj.ecalibModel,cbh);
            % here we update the state variables of the current track
            % according to different types of the kalman filter
            switch trackOption
                case {'force', 'acceleration'}
                    tracks(i).predSpeed = [tracks(i).kalmanFilter.State(2),tracks(i).kalmanFilter.State(5)];
                case 'speed'
                    tracks(i).predSpeed = [tracks(i).kalmanFilter.State(2),tracks(i).kalmanFilter.State(4)];
                otherwise
            end
            
            % here to update bounding boxes, we get a local copy on the
            % location of tracks
            predictedCentroid = tracks(i).predCentroid;
            
            % shift the bounding box so that its center is at
            % the predicted location
            predictedCentroid = int32(predictedCentroid(1:2)) - (int32(bbox(3:4)) / 2);
            
            % reform  the bounding boxes with their new locations
            tracks(i).bbox = double([predictedCentroid(1:2), bbox(3:4)]);
            
            % here to update bounding boxes (on the projected plane), we get a local copy on the
            % location of tracks
            predictedCentroidProj = tracks(i).predCentroidProj;
            
            % shift the bounding box so that its center is at
            % the predicted location
            switch obj.ModPlayerOption
                case 'rectangle'
                    predictedCentroidProj = int32(predictedCentroidProj(1:2)) - (int32(bboxProj(3:4)) / 2);
                    predictedCentroidProj(2) = predictedCentroidProj(2)-bboxProj(4);
                    tracks(i).bboxProj = double([predictedCentroidProj(1:2), bboxProj(3:4)]);
                case 'circle'
                    tracks(i).bboxProj = double([predictedCentroidProj(1:2), bboxProj(3)]);
                otherwise
            end
        end
    end

    function [assignments, unassignedTracks, unassignedDetections] = ...
            detectionToTrackAssignment()
        % in this function we find the assignment between the existing
        % tracks and the observations coming from a certain frame. We do
        % the assignment on the image plane, so the numbers correspond to
        % pixels mostly.
        
        nTracks = length(tracks);
        nDetections = size(centroidsProj, 1);
        % we initialize the assignment matrix and all other related
        % variables first
        if(nTracks == 0 || size(centroidsProj,1) == 0)
            assignments = [];
            unassignedTracks = [];
            unassignedDetections = [1:nDetections];
            return;
        end
        
        % compute the cost of assigning each detection to each track
        % here we want to obtain (x_i - y_j)^2 where x_i is the position of
        % a track and y_j is the position of an observation. for this
        % purpose we first get the x_i^2, then we get the y_j^2, after we
        % get the x_i*y_j. after those it all comes down to calculating
        % x_i^2 - 2*x_i*y_j + y_j^2 using the matrices.
        cent_tracks = cat(1,tracks(:).centroidProj);
        match_mat1 = dot(cent_tracks,cent_tracks,2); % square of tracks x_i^2
        match_mat2 = dot(centroidsProj,centroidsProj,2); % square of detections y_j^2
        match_mat11 = repmat(match_mat1,1,length(match_mat2)); % get x_i^2 into matrix
        match_mat22 = repmat(match_mat2',length(match_mat1),1); % get y_j^2 into matrix
        match_mat12 = cent_tracks*centroidsProj'; % get x_i*y_j into matrix
        
        % we construct the cost matrix here. entry i j refers to the cost
        % of assigning the i'th track to the j'th observation
        cost = match_mat11 - 2*match_mat12 + match_mat22; % get x_i^2 - 2*x_i*y_j + y_j^2 using those matrices
        cost(cost>3500) = inf; % here we do a sanity check on the costs. if the cost of assigning is greater than 2500 meaning
        % if a track and an observation are apart from each other by more than 50 pixels on the image plane,
        % we eliminate the possibilty of an assignment between these two. 2500 is a hand-tuned value.
        
        % the hungarian algorithm for track to observation assignment
        %[matching,~] = Hungarian(cost);
        %matching_inds = find(matching);
        %matched_costs = cost(matching_inds);
        %matching_inds = matching_inds(matched_costs < 2500);
        %[assignments(:,1),assignments(:,2)] = ind2sub(size(matching),matching_inds);
        
        % use of munkres algorithm (a modification of the hungarian
        % algorithm) to do the assignment between tracks and observations
        [matching,~] = munkres(cost);
        matching = matching'; %  get the matching matrix.
        id_tracks = [1:nTracks]';
        assignments = [id_tracks(matching~=0),matching(matching~=0)]; % here we get the assigned tracks
        unassignedTracks = setdiff(id_tracks,assignments(:,1)); % get the tracks that are unassigned
        unassignedDetections = setdiff([1:nDetections],assignments(:,2)); % get the observations that are not assigned
    end

    function updateBlindTracks()
        % this function is used for the case where the kalman filter update
        % does not happen for all the cycles
        numAssignedTracks = size(assignments, 1);
        % the statistical values of all the tracks are updated
        for i = 1:numAssignedTracks % all the tracks that are assigned in the previous step
            trackIdx = assignments(i, 1); % get the track id
            
            if(trackIdx > length(tracks))
                break;
            end
            tracks(trackIdx).age = tracks(trackIdx).age + 1; % update track's age
            tracks(trackIdx).consecutiveInvisibleCount = ...
                tracks(trackIdx).consecutiveInvisibleCount + 1; % update visibility
            tracks(trackIdx).centroid = tracks(trackIdx).kalmanFilter.x'; % update the centroid on the image plane from the kalmen filter
            tracks(trackIdx).centroidProj = tracks(trackIdx).kalmanFilterProj.x'; % update the centroid on the projected plane from the kalmen filter
        end
        for i = 1:length(unassignedTracks) % all the tracks that are not assigned in the previous step
            ind = unassignedTracks(i); % get the track id
            if(ind > length(tracks))
                break;
            end
            tracks(ind).age = tracks(ind).age + 1; % update track's age
            tracks(ind).consecutiveInvisibleCount = ...
                tracks(ind).consecutiveInvisibleCount + 1; % update visibility
            tracks(ind).centroid = tracks(ind).kalmanFilter.x'; % update the centroid on the image plane from the kalmen filter
            tracks(ind).centroidProj = tracks(ind).kalmanFilterProj.x'; % update the centroid on the projected plane from the kalmen filter
        end
    end

    function updateAssignedTracks()
        % here we update then statistical parameters corresponding to
        % assgined tracks along with the kalman filters of those tracks
        % with the observartions associated
        numAssignedTracks = size(assignments, 1);
        for i = 1:numAssignedTracks % for each assignment
            trackIdx = assignments(i, 1); % get the track that is assigned
            detectionIdx = assignments(i, 2); % get the detection id that is assigned
            centroid = centroids(detectionIdx, :); % get the observation correspoding to the position of the detected object, on image plane
            centroidProj = centroidsProj(detectionIdx, :); % get the observation correspoding to the position of the detected object, on projected plane
            bbox = bboxes(detectionIdx, :); % get the observation correspoding to the bounding box of the detected object, on image plane
            bboxProj = bboxesProj(detectionIdx, :); % get the observation correspoding to the bounding box of the detected object, on projected plane
            area = areas(detectionIdx, :); % get the observation correspoding to the area of the detected object, on projected plane
            contour = vertcat(contours{1,detectionIdx}{1,:});
            contour = contour - repmat(centroid,size(contour,1),1);
            contourProj = vertcat(contoursProj{1,detectionIdx}{1,:});
            contourProj = contourProj - repmat(centroidProj,size(contourProj,1),1);
            
            tracks(trackIdx).contour = contour;
            tracks(trackIdx).contourProj = contourProj;
            
            tracks(trackIdx).prevCentroid = tracks(trackIdx).centroid; % get a recording on the previous position of the detected object, on image plane
            tracks(trackIdx).prevCentroidProj = tracks(trackIdx).centroidProj; % get a recording on the previous position of the detected object, on projected plane
            
            % correct the estimate of the object's location
            % using the new detection
            
            % filter out the position observation
            tracks(trackIdx).positionFilter = filterAdd(centroid,tracks(trackIdx).positionFilter);
            pos_t = tracks(trackIdx).positionFilter.value;
            
            % filter out the area observation
            tracks(trackIdx).massFilter = filterAdd(area,tracks(trackIdx).massFilter);
            tracks(trackIdx).area = tracks(trackIdx).massFilter.value;
            
            % fiter out the centroid information
            tracks(trackIdx).centroidFilter = filterAdd(centroid,tracks(trackIdx).centroidFilter);
            tracks(trackIdx).centroid = tracks(trackIdx).centroidFilter.value;
            
            % fiter out the projected centroid information
            tracks(trackIdx).centroidProjFilter = filterAdd(centroidProj,tracks(trackIdx).centroidProjFilter);
            tracks(trackIdx).centroidProj = tracks(trackIdx).centroidProjFilter.value;
            
            % get the observed speed wrt the previous step using a first
            % order differential approximation
            tracks(trackIdx).obsSpeed = tracks(trackIdx).centroid - tracks(trackIdx).prevCentroid;
            
            % update the kalman filter corresponding to the centroid, on
            % the image plane
            tracks(trackIdx).kalmanFilter = tracks(trackIdx).kalmanFilter.Update(centroid',double(tracks(trackIdx).area));
            tracks(trackIdx).centroid = tracks(trackIdx).kalmanFilter.x';
            
            % update the kalman filter corresponding to the centroid, on
            % the projected
            tracks(trackIdx).kalmanFilterProj = tracks(trackIdx).kalmanFilterProj.Update(centroidProj',double(tracks(trackIdx).area));
            tracks(trackIdx).centroidProj = tracks(trackIdx).kalmanFilterProj.x';
            % for some certaion kind of the kalman filter depending on the
            % model based on the force applied on the mass this update is
            % used
            %tracks(trackIdx).kalmanFilter = tracks(trackIdx).kalmanFilter.Update(pos_t',double(tracks(trackIdx).area));
            
            %             if(strcmp(trackOption,'areas'))
            %                 correct(tracks(trackIdx).kalmanFilter, [tracks(trackIdx).centroid,tracks(trackIdx).area]);
            %             else
            %                 correct(tracks(trackIdx).kalmanFilter, tracks(trackIdx).centroid);
            %             end
            
            % getting the speed estimate from the kalman depending on the type of kalman used...
            switch trackOption
                case {'force','acceleration'}
                    speed = [tracks(trackIdx).kalmanFilter.State(2), tracks(trackIdx).kalmanFilter.State(5)];
                case 'speed'
                    speed = [tracks(trackIdx).kalmanFilter.State(2), tracks(trackIdx).kalmanFilter.State(4)];
                otherwise
            end
            
            % filter out the speed information
            tracks(trackIdx).speedFilter = filterAdd(speed, tracks(trackIdx).speedFilter);
            speed = tracks(trackIdx).speedFilter.value;% * camParams.conversionFactor;
            
            % for some reason, in case needed, get the x and y speeds
            tracks(trackIdx).xspeed = speed(2);
            tracks(trackIdx).yspeed = speed(1);
            
            % replace predicted bounding box with detected
            % bounding box
            tracks(trackIdx).bbox = bbox;
            tracks(trackIdx).bboxProj = bboxProj;
            
            % update track's age
            tracks(trackIdx).age = tracks(trackIdx).age + 1;
            
            % update visibility
            tracks(trackIdx).totalVisibleCount = ...
                tracks(trackIdx).totalVisibleCount + 1;
            
            % since we know that the track was visible for the previous
            % step, we set this variable to zero, as the name also suggests
            tracks(trackIdx).consecutiveInvisibleCount = 0;
        end
    end

    function updateUnassignedTracks()
        % here we update the statistical parameters corresponding to
        % unassigned tracks
        for i = 1:length(unassignedTracks) % for each unassigned track we loop
            ind = unassignedTracks(i); % get the id of the track
            tracks(ind).age = tracks(ind).age + 1; % increment the age by one
            tracks(ind).consecutiveInvisibleCount = ...
                tracks(ind).consecutiveInvisibleCount + 1; % since it is not observed, increase the invisible count
        end
    end

    function deleteLostTracks()
        if isempty(tracks)
            return;
        end
        
        invisibleForTooLong = 4;
        ageThreshold = 3;
        
        % compute the fraction of the track's age for which it was visible
        ages = [tracks(:).age];
        totalVisibleCounts = [tracks(:).totalVisibleCount];
        visibility = totalVisibleCounts ./ ages;
        
        % find the indices of 'lost' tracks
        lostInds = (ages < ageThreshold & visibility < 0.1) | ...
            [tracks(:).consecutiveInvisibleCount] >= invisibleForTooLong;
        
        % delete lost tracks
        tracks = tracks(~lostInds);
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

    function createNewTracks()
        % here we create new tracks for the obbservations that are not
        % assigned to any existing tracks
        
        % first we filter out the necessary information corresponding to
        % unassigned observations. kind of information should be obvious
        % from looking at the name
        centroidsU = centroids(unassignedDetections, :);
        centroidsProjU = centroidsProj(unassignedDetections, :);
        bboxesU = bboxes(unassignedDetections, :);
        bboxesProjU = bboxesProj(unassignedDetections, :);
        areasU = areas(unassignedDetections, :);
        contoursU = contours(unassignedDetections);
        contoursProjU = contoursProj(unassignedDetections);
        % for each unassigned observation we loop
        for i = 1:size(centroidsU, 1)
            
            % we get the information correponding to the ith observation.
            % kind of inforamtion should be obvious looking at the name
            centroid = centroidsU(i,:);
            centroidProj = centroidsProjU(i,:);
            bbox = bboxesU(i, :);
            bboxProj = bboxesProjU(i, :);
            contour = vertcat(contoursU{1,i}{1,:});
            contour = contour - repmat(centroid,size(contour,1),1);
            contourProj = vertcat(contoursProjU{1,i}{1,:});
            contourProj = contourProj - repmat(centroidProj,size(contourProj,1),1);
            
            area = areasU(i);
            
            % create the kalman filter with the initial parameters. Order
            % and the correspondance of the inputs are given by the below
            % signature of the initialization function. This filter is used
            % for tracking purposes on the image plane.
            % SimpleKF(init_x,sigma_p,sigma_q,sigma_r,type,m,time)
            kalmanFilter = SimpleKF(centroid,0.6,[0.01,0.01],1e1,trackOption,double(area),1);
            
            % to track on the projected plane, we initialize another kalman
            % filter. The order of the input are the same as the previous
            % one.
            kalmanFilterProj = SimpleKF(centroidProj,0.6,[0.01,0.01],1e1,trackOption,double(area),1);
            
            % create a new track
            newTrack = struct(...
                'id', nextId, ... % id of the track
                'bbox', bbox, ... % the bounding box corresponding to object, on the image axis
                'bboxProj', bboxProj, ... % the bounding box corresponding to object, projected onto real-world axis
                'kalmanFilter', kalmanFilter, ... % kalman filter corresponding to this object
                'kalmanFilterProj', kalmanFilterProj, ... % kalman filter corresponding to this object, projected onto real-world axis
                'age', 1, ... % how old the track is
                'totalVisibleCount', 1, ... % for how many cycles do we get an observation on this track
                'consecutiveInvisibleCount', 0, ... % if we cannot observe it, for how many previous steps we could not
                'xspeed', 0, ... % a variable for keeping one directional speed in image plane
                'yspeed', 0, ... % a variable for keeping one directional speed in image plane
                'predSpeed', [0,0], ... % a variable keeping the 2D predicted speed in image plane
                'obsSpeed', [0,0], ... % a variable keeping the 2D observed speed in image plane
                'area', double(area), ... % a variable keeping the area of the tracked object
                'contour', contour, ... % a variable keeping the relevant contour points on the image plane
                'contourProj', contourProj, ... % a variable keeping the relevant contour points on the real world coordinates
                'speedFilter', filterInit(10,2), ... % a filter for the observed speed on the image plane
                'positionFilter', filterInit(10,2), ... % a filter for the position on the image plane
                'centroidFilter', filterInit(1,2), ... % a filter for the centroid on the image plane
                'centroidProjFilter', filterInit(1,2), ... % a filter for the centroid on the real plane
                'massFilter', filterInit(10,1), ... % a filter on the predicted mass (not used)
                'centroid', centroid, ... % the representative position of the center of the object on the image plane
                'centroidProj', centroidProj, ... % the representative position of the center of the object on the projected plane
                'predCentroid', centroid, ... % prediction on the representative position of the center of the object on the image plane
                'predCentroidProj', centroidProj, ... % prediction on the representative position of the center of the object on the projected plane
                'prevCentroid', centroid, ... % previous centroid position on the image plane
                'prevCentroidProj', centroidProj); % previous centroid position on the projected plane
            
            
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
        if(obj.track)
            minVisibleCount = 1;
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
                    realKals = cat(1,reliableTracks.kalmanFilterProj);
                    %!!!!! pay attention !!!!!
                    realProjs = cat(1,realKals.x);
                    ids = cat(1,reliableTracks.id);
                    % get ids
                    %                 ids = int32([reliableTracks(:).id]);
                    %                 xspeeds = [reliableTracks(:).xspeed];
                    %                 yspeeds = [reliableTracks(:).yspeed];
                    
                    % create labels for objects indicating the ones for
                    % which we display the predicted rather than the actual
                    % location
                    %                 labels = cellstr(int2str(ids'));
                    %                 xlspeeds = cellstr(num2str(xspeeds',2));
                    %                 ylspeeds = cellstr(num2str(yspeeds',2));
                    
                    %predictedTrackInds = ...
                    %    [reliableTracks(:).consecutiveInvisibleCount] > 0;
                    %isPredicted = cell(size(labels));
                    %isPredicted(predictedTrackInds) = {' predicted'};
                    %space = {', '};
                    %                 labels = strcat(labels, ', x: ', xlspeeds, ', y: ', ylspeeds);
                    
                    predictedTrackInds = ...
                        [reliableTracks(:).consecutiveInvisibleCount] > 0;
                    
                    newTrackInds = ...
                        [reliableTracks(:).consecutiveInvisibleCount] <= 0 & [reliableTracks(:).age] <= 5;
                    trackedTrackInds = ...
                        [reliableTracks(:).consecutiveInvisibleCount] <= 0 & [reliableTracks(:).age] > 5;
                    
                    if (sum(predictedTrackInds)>0)
                        bboxColor = 'r';
                        % draw on the mask
                        mask = objectAnnotation(mask, 'Rectangle', bboxes(predictedTrackInds,:),bboxColor,3);
                        %mask = objectAnnotation(mask, 'Text', [ids(predictedTrackInds,:),bboxes(predictedTrackInds,1:2)],bboxColor,3);
                    end
                    if (sum(trackedTrackInds)>0)
                        bboxColor = 'y';
                        % draw on the mask
                        mask = objectAnnotation(mask, 'Rectangle', bboxes(trackedTrackInds,:),bboxColor,3);
                        %mask = objectAnnotation(mask, 'Text', [ids(trackedTrackInds,:),bboxes(trackedTrackInds,1:2)],bboxColor,3);
                    end
                    if (sum(newTrackInds)>0)
                        bboxColor = 'g';
                        % draw on the mask
                        mask = objectAnnotation(mask, 'Rectangle', bboxes(newTrackInds,:),bboxColor,3);
                        %mask = objectAnnotation(mask, 'Text', [ids(newTrackInds,:),bboxes(newTrackInds,1:2)],bboxColor,3);
                    end
                    
                    % draw on the mask
                    %mask = objectAnnotation(mask, 'Rectangle', bboxes,bboxColor,4);
                    
                    % draw on the frame
                    frame = cv.drawContours(frame,contours,'Color',[255 255 0],'Thickness',3);
                    predFrame = cv.drawContours(predFrame,contours,'Color',[0 255 0],'Thickness',3);
                end
            end
            allPredContours = [predCloudContours,predShadowContours];
            colors = [repmat('b',size(predCloudContours)),repmat('g',size(predShadowContours))];
            % draw on the mask
            
            %mask = objectAnnotation(mask, 'Rectangle', bboxes);
            
            % draw on the frame
            %frame = cv.drawContours(frame,contours,'Color',[255 255 0]);
            
            % display the mask and the frame
            obj.maskPlayer.Step(mask);
            obj.videoPlayer.Step(frame);
            obj.predPlayer.Step(predFrame);
            obj.mapPlayer = obj.mapPlayer.Step(allPredContours,colors);
            if(obj.gifMode)
                frame1 = obj.videoPlayer.GetFrame();
                frame2 = obj.maskPlayer.GetFrame();
                frame3 = obj.predPlayer.GetFrame();
                pause(0.02);
                frame4 = obj.mapPlayer.GetFrame();
                ff = [frame1,frame2;frame3,frame4];
                if(mod(obj.frameCount-obj.frameRate,50*obj.frameRate)<=obj.frameRate)
                    [ffind,cm] = rgb2ind(ff,256);
                    imwrite(ffind,cm,strcat(obj.frameRecordImageName,'.gif'),'gif', 'Loopcount',inf);
                else
                    [ffind,cm] = rgb2ind(ff,256);
                    imwrite(ffind,cm,strcat(obj.frameRecordImageName,'.gif'),'gif', 'WriteMode','append');
                end
            end
        else
            % if we are dumping to a webseite.
            textColor = 'g';
            irradianceValues = obj.irradianceValues(end);
            value_positions = [repmat(70,1,length(irradianceValues));70:35:70+((length(irradianceValues)-1)*35)];
            frame = objectAnnotation(frame, 'Text', [irradianceValues, value_positions'], textColor,3);
            obj.videoPlayer.Step(frame,obj.frameDateStr);
            obj.plotPlayerGHI.Step([obj.irradianceTimes,obj.irradianceValues],strcat('Time Relative to ', obj.frameDateStr, ' (Min.)') ,'Irradiance (W/m2)','Irradiance Value Changes');
            obj.plotPlayerPV.Step([obj.powerTimes,obj.powerValues],strcat('Time Relative to ', obj.frameDateStr, ' (Min.)') ,'Power (W)','Power Production Value Changes');
            obj.plotPlayerTemp.Step([obj.tempTimes,obj.tempValues],strcat('Time Relative to ', obj.frameDateStr, ' (Min.)') ,'Temperature (degrees)','Temperature Value Changes');
            
            if(obj.gifMode)
                frame1 = obj.videoPlayer.GetFrame();
                frame2 = obj.plotPlayerGHI.GetFrame();
                frame3 = obj.plotPlayerPV.GetFrame();
                frame4 = obj.plotPlayerTemp.GetFrame();
                ff = [frame1,frame2;frame3,frame4];
                %imwrite(ff,strcat(obj.frameRecordName,'.jpg'));
                if(mod(obj.frameCount-1,20)==1)
                    [ffind,cm] = rgb2ind(ff,256);
                    imwrite(ffind,cm,strcat(obj.frameRecordImageName,'.gif'),'gif', 'Loopcount',inf);
                    fid = fopen(strcat(obj.frameRecordPageName,'.html'),'w');
                    fprintf(fid, '%s\n', strcat('<img src="', strcat(obj.frameRecordImagePageName,'.gif'), '" alt="the file here">'));
                    fclose(fid);
                else
                    [ffind,cm] = rgb2ind(ff,256);
                    imwrite(ffind,cm,strcat(obj.frameRecordImageName,'.gif'),'gif', 'WriteMode','append');
                end
            else
                if(mod(obj.frameCount-1,20)==1)
                    frame1 = obj.videoPlayer.GetFrame();
                    frame2 = obj.plotPlayerGHI.GetFrame();
                    frame3 = obj.plotPlayerPV.GetFrame();
                    frame4 = obj.plotPlayerTemp.GetFrame();
                    ff = [frame1,frame2;frame3,frame4];
                    imwrite(ff,strcat(obj.frameRecordImageName,'.jpg'));
                    fid = fopen(strcat(obj.frameRecordPageName,'.html'));
                    fprintf(fid, '%s\n', strcat('<img src="', strcat(obj.frameRecordImagePageName,'.jpg'), '" alt="the file here">'));
                    fclose(fid);
                end
                
            end
            
        end
    end

    function copyToWebServer()
    end
%% Generic filter implementation.
% c : length of the filter
% d : dimensions of the filter
% functions for initializing a filter, we can build a separate
% class out of this part...
% B. Zeydan, 27. Feb. 2013
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
        %filter.value = value';
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
% B. Zeydan, 27. Feb. 2013
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
        
        cbh = 300;
        
        timeConversionFactor = 1;
        lengthConversionFactor = cbh / sqrt(refsensor^2  + ((xefsensor/2)/tan(fov/2))^2);
        conversionFactor = 1;%rpixel * lengthConversionFactor / timeConversionFactor;
        
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