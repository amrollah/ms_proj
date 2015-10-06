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

function [meanCentroids,meanCentroidsProj,centroids,centroidsProj]=multiObjectTrackingCV_CamCalibAssistedLabel(obj,tracksRecorded,track_id)
% create system objects used for 2reading video, detecting moving objects,
% and displaying the results

meanCentroids = [];
meanCentroidsProj = [];
tracks = initializeTracks(); % create an empty array of tracks
contours = [];
centroids = [];
centroidsProj = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the variables that determine the options and parameters applied on the
% tracking algorithm.

cbh = 900; % the assumed cloud base height parameter that determines the projection performance of the algorithm.
ModPlayerOption = 'map'; % the representation of clouds on the rea-sized model is done using this opttion. rectangle is another possibility.
maskMethod = 'body'; % identification type of clouds is done by this variable. either they are identified by their body or the contours (edges)
trackOption = 'speed'; % the option that is applied for the kalman filter, for other possibilities please look inside the object implementation.

sunEliminationCycle = 100; % every this many frames, we calculate the sun's position once again.
sunPosition = [-1,-1]; % we initialize the position of sun, to be used throughout the code
deg2radf = pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% in case we go with the validation mode and we are recording the
% performance of the algorithm, the pstate variables of all the steps are
% recorded in these variables
obj = setupSystemObjects(obj);

if(isempty(obj.reader))
    disp('Camera does not respond. Terminating...');
    return;
end

%writerObj = VideoWriter(['Validation_Sequence_',num2str(track_id),'.avi']);
%writerObj.FrameRate = 1;
%open(writerObj);
% This is the main loop of the program
itr = 1;
itr2 = 1;
label = 1;
while obj.frameCount < obj.finalFrame
    if(label)
        if(mod(itr,50)==1)
            itr
        end
        tracks = cell2mat(tracksRecorded(itr));
        ids = [tracks(:).id];
        track = tracks(ids==track_id);

        if(sum(ids == track_id))
            frame = readFrame();
            contours = labelTheImage(track);
            if(sum(track.centroidObs)~=-2)
                meanCentroids(itr2,:) = mean(track.contour+1 + repmat(track.centroidObs,size(track.contour,1),1),1); 
                meanCentroidsProj(itr2,:) = mean(track.contourProj + repmat(track.centroidProjObs,size(track.contourProj,1),1), 1);
                centroids(itr2,:) = track.centroidObs;
                centroidsProj(itr2,:) = track.centroidProjObs;
                itr2 = itr2+1;
            end
            displayLabelingResults();
            
            %writeVideo(writerObj,frame);
        else
            obj.frameCount = obj.frameCount+obj.frameRate;
        end
        
    else
        frame = readFrame();
        tracks = cell2mat(tracksRecorded(itr));
        ids = [tracks(:).id];
        contours = getContours(tracks);
        for ir = 1:length(ids)
            frame = cv.putText(frame,num2str(ids(ir)),tracks(ir).centroidObs,'Color',[1,1,0],'Thickness',4);
        end
        displayLabelingResults;
    end
    itr = itr+1;
    %evaluateTrackingData(obj.timeHorizon);
    %predictNewLocationsOfTracksIntoFuture(obj.timeHorizon);
end
%close(writerObj);

    function obj = setupSystemObjects(obj)
        % Initialize Video I/O
        % Create objects for reading a video from a file, drawing the tracked
        % objects in each frame, and playing the video.
        obj.frameCount = 1;
        obj.videoPlayer = VidPlayer([20, 500, 700, 400]); % for the real image
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
            'centroidObs', {}, ... % the representative observed position of the center of the object on the image plane
            'centroidProj', {}, ... % the representative position of the center of the object on the projected plane
            'centroidProjObs', {}, ... % the representative observed position of the center of the object on the projected plane
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

    function [contours, contoursProj, areas, areasProj, centroids, centroidsProj, bboxes, bboxesProj, mask] = detectObjects(frame)
        
        % detect objects into a binary image
        mask = getBinaryImage(frame,maskMethod);
        mask = medfilt2(mask,[10 10]);
        % apply morphological operations to remove noise and fill in holes
        %mask = imopen(mask, strel('disk', 3));
        %mask = imclose(mask, strel('square', 30));
        %mask = imfill(mask, 'holes');
        %contours = cv.findContours(mask,'Mode','List'); % contour extraction using openCV method
        contours = cv.findContours(mask,'Mode','External'); % contour extraction using openCV method
        % areas of objects are approximately computed here using several methods
        %areas = cell2mat(cellfun(@(x) size(cell2mat(x(:)),1),contours,'UniformOutput',false))';
        %areas = cell2mat(cellfun(@(x) prod(max(cell2mat(x(:)),[],1)-min(cell2mat(x(:)),[],1)) ,contours,'UniformOutput',false))';
        areas = cell2mat(cellfun(@(x) cv.contourArea(x(:)) ,contours,'UniformOutput',false))';
        
        % elminate the objects with contour smaller than 350.
        consider = (areas > 100);% & (areas < 100500); % if the area of an object is less than 350 pixels, we discard that
        
        % take only the objects that we want to consider
        contours = contours(consider);
        areas = areas(consider);
        %+ones(size(cell2mat(x(:))))
        contoursProj = cellfun(@(x) num2cell(getImage2Real(cell2mat(x(:))+1,obj.calibModel,obj.ecalibModel,obj.cbh),2)',contours,'UniformOutput',false);
        areasProj = cell2mat(cellfun(@(x) cv.contourArea(x(:)) ,contoursProj,'UniformOutput',false))';
        % different methods for determining the point that represent the
        % cloud objets first method uses the mean of the contours the
        % second one uses the mid point of the bounding rectangle
        %centroids = cellfun(@(x) mean(cell2mat(x(:)),1),contours,'UniformOutput',false);
        centroids = cellfun(@(x) (min(cell2mat(x(:))+1,[],1)+max(cell2mat(x(:))+1,[],1))/2,contours,'UniformOutput',false);
        centroidsProj = cellfun(@(x) (min(cell2mat(x(:)),[],1)+max(cell2mat(x(:)),[],1))/2,contoursProj,'UniformOutput',false);
        
        bboxes = cellfun(@(x) [min(cell2mat(x(:))+1,[],1),max(cell2mat(x(:))+1,[],1)-min(cell2mat(x(:))+1,[],1)] ,contours,'UniformOutput',false);
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
        frame = cv.circle(frame,[x,y],90,'Color',[0,0,0],'Thickness',-1);
    end

    function getShadowDisplacement(model3D)
        DN = datenum(obj.frameDateVec(1:6)');
        Time = pvl_maketimestruct(DN, model3D.UTC);
        [SunAz, ~, ApparentSunEl, ~]=pvl_ephemeris(Time, model3D.Location);
        zenith = 90-ApparentSunEl;
        azimuth = SunAz;
        [x,y,~] = sph2cart((90-azimuth)*deg2radf,(90-zenith)*deg2radf,obj.cbh/cos(zenith*deg2radf));
        obj.shadowDisp = [-x,-y];
    end

    function tracksPred = evaluateTrackingData(timeSteps)
        %tracksPred = [tracks,tracksLost];
        tracksPred = cell2mat(tracks);
        for i = 1:length(tracksPred)
            for predNum = 1:timeSteps
                if(length(tracksRecorded)<itr+predNum)
                    break;
                end
                trackObserved = tracksRecorded{itr+predNum}([tracksRecorded{itr+predNum}(:).id]==tracksPred(i).id);
                if(isempty(trackObserved))
                    break;
                end
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
                
                predictionErrors(itr,tracksPred(i).id,predNum).centroidObs = trackObserved.centroidObs;
                predictionErrors(itr,tracksPred(i).id,predNum).centroidPred = tracksPred(i).predCentroid;
                predictionErrors(itr,tracksPred(i).id,predNum).errCentroid = trackObserved.centroidObs - tracksPred(i).predCentroid;
                predictionErrors(itr,tracksPred(i).id,predNum).errCentroidN = norm(trackObserved.centroidObs - tracksPred(i).predCentroid);
                predictionErrors(itr,tracksPred(i).id,predNum).centroidProjObs = trackObserved.centroidProjObs;
                predictionErrors(itr,tracksPred(i).id,predNum).centroidProjPred = tracksPred(i).predCentroidProj;
                predictionErrors(itr,tracksPred(i).id,predNum).errCentroidProj = trackObserved.centroidProjObs - tracksPred(i).predCentroidProj;
                predictionErrors(itr,tracksPred(i).id,predNum).errCentroidProjN = norm(trackObserved.centroidProjObs - tracksPred(i).predCentroidProj);
                if(sum(trackObserved.centroidObs)==-2)
                    %predictionErrors(itr,tracksPred(i).id,predNum).centroidObs = [];
                    predictionErrors(itr,tracksPred(i).id,predNum).errCentroid = [];
                    predictionErrors(itr,tracksPred(i).id,predNum).errCentroidN = [];
                end
                if(sum(trackObserved.centroidProjObs)==-2)
                    %predictionErrors(itr,tracksPred(i).id,predNum).centroidProjObs = [];
                    predictionErrors(itr,tracksPred(i).id,predNum).errCentroidProj = [];
                    predictionErrors(itr,tracksPred(i).id,predNum).errCentroidProjN = [];
                end
            end
        end
    end

    function tracksPred = predictNewLocationsOfTracksIntoFuture(timeSteps)
        %tracksPred = [tracks,tracksLost];
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

    function labelContour = labelTheImage(track)
        labelContour = [];
        if(norm(track.centroidObs-[-1,-1]))
            labelContour{1,1} = num2cell(getReal2Image(track.contourProj+repmat(track.centroidProjObs,size(track.contourProj,1),1),obj.calibModel,obj.ecalibModel,obj.cbh)-1,2)';
        end
    end
    
    function predContours = getContours(tracksPred)
        predContours=[];
        iterator = 1;
        for i = 1:length(tracksPred)
            if(tracksPred(i).consecutiveInvisibleCount<1)
                predContours{1,iterator} = num2cell(getReal2Image(tracksPred(i).contourProj+repmat(tracksPred(i).centroidProjObs,size(tracksPred(i).contourProj,1),1),obj.calibModel,obj.ecalibModel,obj.cbh)-1,2)';
                iterator = iterator+1;
            end
        end
    end

    function predFrame = getPredictedImage(tracksPred)
        predFrame = zeros(size(frame));
        predContours=[];
        iterator = 1;
        for i = 1:length(tracksPred)
            if(tracksPred(i).consecutiveInvisibleCount<1)
                predContours{1,iterator} = num2cell(getReal2Image(tracksPred(i).contourProj+repmat(tracksPred(i).predCentroidProj,size(tracksPred(i).contourProj,1),1),obj.calibModel,obj.ecalibModel,obj.cbh)-1,2)';
                iterator = iterator+1;
            end
        end
        if(length(tracksPred)~=0)
            predFrame = cv.drawContours(predFrame,predContours,'Thickness',-1,'LineType',4);
        end
    end

    function [occluded] = getOcclusion(tracksPred,model3D)
        occluded = 0;
        for i = 1:length(tracksPred)
            shadowContProj  = tracksPred(i).contourProj+repmat(tracksPred(i).predCentroidProj+obj.shadowDisp,size(tracksPred(i).contourProj,1),1);
            [in,on] = inpolygon(model3D.Sen_xs,model3D.Sen_ys,shadowContProj(:,1),shadowContProj(:,2));
            occluded = occluded | in | on;
        end
    end

    function [predCloudContours,predShadowContours] = getPredictedScenario(tracksPred)
        predCloudContours=[];
        predShadowContours=[];
        for i = 1:length(tracksPred)
            predCloudContours{1,i} = num2cell([tracksPred(i).contourProj+repmat(tracksPred(i).predCentroidProj,size(tracksPred(i).contourProj,1),1),obj.cbh*ones(size(tracksPred(i).contourProj,1),1)],2)';
            predShadowContours{1,i} = num2cell([tracksPred(i).contourProj+repmat(tracksPred(i).predCentroidProj+obj.shadowDisp,size(tracksPred(i).contourProj,1),1),zeros(size(tracksPred(i).contourProj,1),1)],2)';
            %predShadowContoursP{1,i} = num2cell([tracksPred(i).contourProj+repmat(tracksPred(i).predCentroidProj+obj.shadowDisp+[model3D.Cam_xc,model3D.Cam_yc],size(tracksPred(i).contourProj,1),1)],2)';
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

    function [assignments, unassignedTracks, unassignedDetections] = detectionToTrackAssignment()
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
        cent_tracks = cat(1,tracks(:).predCentroidProj);
        match_mat1 = dot(cent_tracks,cent_tracks,2); % square of tracks x_i^2
        match_mat2 = dot(centroidsProj,centroidsProj,2); % square of detections y_j^2
        match_mat11 = repmat(match_mat1,1,length(match_mat2)); % get x_i^2 into matrix
        match_mat22 = repmat(match_mat2',length(match_mat1),1); % get y_j^2 into matrix
        match_mat12 = cent_tracks*centroidsProj'; % get x_i*y_j into matrix
        
        % we construct the cost matrix here. entry i j refers to the cost
        % of assigning the i'th track to the j'th observation
        cost = match_mat11 - 2*match_mat12 + match_mat22; % get x_i^2 - 2*x_i*y_j + y_j^2 using those matrices
        cost_th = (obj.cbh^2/300^2)*1e4;
        cost(cost>cost_th) = inf; % here we do a sanity check on the costs. if the cost of assigning is greater than 2500 meaning
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

    function [assignments, unassignedTracks, unassignedDetections] = detectionToTrackAssignmentGaussian()
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
        
        
        for i = 1:nTracks
            P = tracks(i).kalmanFilterProj.P;
            P_s = [P(1,1),P(1,3),0;P(3,1),P(3,3),0;0,0,1e6];
            cost(i,:) = 1-mvnpdf([centroidsProj,areas],[tracks(i).predCentroidProj,tracks(i).area],P_s)';
        end
        cost(cost==1)=inf;
        %cost(cost>3500) = inf; % here we do a sanity check on the costs. if the cost of assigning is greater than 2500 meaning
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
            
            tracks(trackIdx).centroidObs = centroid;
            tracks(trackIdx).centroidProjObs = centroidProj;
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
            tracks(trackIdx).kalmanFilter = tracks(trackIdx).kalmanFilter.Update(tracks(trackIdx).centroid',double(tracks(trackIdx).area));
            tracks(trackIdx).centroid = tracks(trackIdx).kalmanFilter.x';
            
            % update the kalman filter corresponding to the centroid, on
            % the projected
            tracks(trackIdx).kalmanFilterProj = tracks(trackIdx).kalmanFilterProj.Update(tracks(trackIdx).centroidProj',double(tracks(trackIdx).area));
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
            tracks(ind).centroid = tracks(ind).predCentroid;
            tracks(ind).centroidObs = [-1,-1];
            tracks(ind).centroidProj = tracks(ind).predCentroidProj;
            tracks(ind).centroidProjObs = [-1,-1];
            tracks(ind).age = tracks(ind).age + 1; % increment the age by one
            tracks(ind).consecutiveInvisibleCount = ...
                tracks(ind).consecutiveInvisibleCount + 1; % since it is not observed, increase the invisible count
        end
    end

    function deleteLostTracks()
        if isempty(tracks)
            return;
        end
        
        invisibleForTooLong = 5;
        ageThreshold = 3;
        
        % compute the fraction of the track's age for which it was visible
        ages = [tracks(:).age];
        totalVisibleCounts = [tracks(:).totalVisibleCount];
        visibility = totalVisibleCounts ./ ages;
        
        % find the indices of 'lost' tracks
        lostInds = (ages < ageThreshold & visibility < 0.1) | ...
            [tracks(:).consecutiveInvisibleCount] >= invisibleForTooLong;
        
        % delete lost tracks
        %tracksLost = [tracksLost,tracks(lostInds)];
        tracks = tracks(~lostInds);
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
        areasProjU = areasProj(unassignedDetections, :);
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
            areaProj = areasProjU(i);
            % create the kalman filter with the initial parameters. Order
            % and the correspondance of the inputs are given by the below
            % signature of the initialization function. This filter is used
            % for tracking purposes on the image plane.
            % SimpleKF(init_x,sigma_p,sigma_q,sigma_r,type,m,time)
            kalmanFilter = SimpleKF(centroid,0.6,[0.01,0.01],1e1,trackOption,double(area),1);
            
            % to track on the projected plane, we initialize another kalman
            % filter. The order of the input are the same as the previous
            % one.
            kalmanFilterProj = SimpleKF(centroidProj,(obj.cbh^2/300^2)*1e2,(obj.cbh^2/300^2)*[9e4,1e3],(obj.cbh^2/300^2)*2e4,trackOption,double(areaProj),1);
            %kalmanFilterProj = SimpleKF(centroidProj,1e2,[9e4,1e3],2e4,trackOption,double(areaProj),1);
            %kalmanFilterProj = SimpleKF(centroidProj,0.6,[0.01,0.01],1e1,trackOption,double(areaProj),1);
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
                'centroidObs', centroid, ... % the representative observed position of the center of the object on the image plane
                'centroidProj', centroidProj, ... % the representative position of the center of the object on the projected plane
                'centroidProjObs', centroidProj, ... % the representative observed position of the center of the object on the projected plane
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

    function displayLabelingResults()
        
        if(~isempty(contours))
            frame = cv.drawContours(frame,contours,'Color',[255 255 0],'Thickness',3);
        end
        obj.videoPlayer.Step(frame);
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
        
        cbh_l = 500;
        
        timeConversionFactor = 1;
        lengthConversionFactor = cbh_l / sqrt(refsensor^2  + ((xefsensor/2)/tan(fov/2))^2);
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
            'cbh', cbh_l, ...
            'xfov', fov);
    end
% end of camera parameter and cbh dependent part...

end