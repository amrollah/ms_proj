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
classdef cloudTrackingCV
    properties (GetAccess=private)
        frame
        predFrame %mask with the predicted future cloud estimates
        predContours %contours of the predicted future cloud estimates. Not implemented here, check multiObjectTrackingCV_CamCalib.m
        mask % the element that keeps the binary image that contains clouds only
        tracks % array of tracks
        tracksPred
        perfs
        nextId % ID of the next track
        centroids % the array that keeps the centroid positions of tracks on the image plane
        centroidsProj % the array that keeps the centroid positions of tracks on the projected (real) plane
        bboxes % the array that keeps the bounding boxes of the detected cloud trails in the image plane
        bboxesProj % the array that keeps the bounding boxes of the detected cloud trails in the projected (real) plane
        contours
        contoursProj
        areas % the array that keeps the areas of bounding boxes in the image plane
        areasProj % the array that keeps the areas of bounding boxes in the projected (real) plane
        predCloudContours
        predShadowContours
        trackOption % the option that is applied for the kalman filter, for other possibilities please look inside the object implementation.
        deg2radf
        shadowDisp
        shadowDisp_d
        live
        vid_r
        vid_n
        vidObj
        reader
        vid_ims
        im_list
        finalFrame
        gifMode % if we want to export minutely recordings as a gif to the web page or not.
        frameCount % the order of frames that are processed.
        frameName % the file name of the frame that is processed the last.
        calibModel % the internal calibration model that is used by the algorithm
        ecalibModel % the external calibration model that is used by the algorithm
        cbh % the assumed cloud base height that is used.
        ModPlayerOption % the way clouds represented in the model player (can be cycle or rectangle)
        predict % whether we are in the tracking mode or display mode. can be 0 or 1
        timeSeriesGHI % the timeseries object corresponding to the GHI measurements
        timeSeriesGHI_cs
        timeSeriesPV % the timeseries object corresponding to the PV measurements
        timeSeriesTemp % the timeseries object corresponding to the temperature measurements
        timeSeriesOcc % the timeseries object corresponding to the occlusion measurements
        timeSeriesCBH % the timeseries object corresponding to the CBH measurements
        irradianceValue % the variable that keeps the irradiance value corresponding to the current frame
        imageMask
        frameRate
        timeHorizon
        irradianceValues
        irradianceTimes
        irradianceValues_cs
        irradianceTimes_cs
        occlusionPoints
        occlusionTimes
        occlusionTimeNums
        occlusionValues
        occlusionQualities
        occValues
        occTimes
        powerValues
        powerTimes
        occState
        clearnessState
        attCs
        attOcc
        irrCs
        irrMs
        irradianceValues_est
        tempValues
        tempTimes
        shadowGrx
        shadowGry
        videoPlayer
        maskPlayer
        modPlayer
        predPlayer
        mapPlayer
        plotPlayerTemp
        plotPlayerPV
        plotPlayerGHI
        frameRecordImageName
        frameRecordImagePageName
        frameRecordPageName
        frameNameStr
        frameDateStr
        frameDateVec
        sunEliminationCycle
        sunPosition
        sunPatch
        sunPattern
        X
        Y
        maskMethod
        recording_tracks
        model3D
        assignments
        unassignedTracks
        unassignedDetections
        unassignedDetectionsToTracks
        unassignedTracksToDetections
    end
    methods
        function obj = cloudTrackingCV(vid_n,vid_r,image_folder,live,verbose,calibModel,ecalibModel,model3D,...
                imageMask,modelImage,cbh,ModPlayerOption,predict,timeSeriesGHI,timeSeriesGHI_cs,...
                timeSeriesPV,timeSeriesTemp,timeSeriesOcc,timeSeriesCBH,maskMethod,trackOption, sunPattern)
            % Initialize Video I/O
            % Create objects for reading a video from a file, drawing the tracked
            % objects in each frame, and playing the video.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % initialization of the key variables that are used as part of a model
            obj.tracks = cloudTrackingCV.initializeTracks(); % create an empty array of tracks
            obj.perfs = cloudTrackingCV.initializePerfs();
            %tracksLost = initializeTracks();
            obj.nextId = 1; % ID of the next track
            obj.centroids = []; % the array that keeps the centroid positions of tracks on the image plane
            obj.centroidsProj = []; % the array that keeps the centroid positions of tracks on the projected (real) plane
            obj.bboxes = []; % the array that keeps the bounding boxes of the detected cloud trails in the image plane
            obj.bboxesProj = []; % the array that keeps the bounding boxes of the detected cloud trails in the projected (real) plane
            obj.areas = []; % the array that keeps the areas of bounding boxes
            obj.areasProj = []; % the array that keeps the areas of bounding boxes
            obj.mask = []; % the element that keeps the binary image that contains clouds only
            obj.predCloudContours = [];
            obj.predShadowContours = [];
            obj.recording_tracks = [];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % the variables that determine the options and parameters applied on the
            % tracking algorithm.
            obj.maskMethod = maskMethod; % identification type of clouds is done by this variable. either they are identified by their body or the contours (edges)
            obj.trackOption = trackOption; % the option that is applied for the kalman filter, for other possibilities please look inside the object implementation.
            obj.sunEliminationCycle = 100; % every this many frames, we calculate the sun's position once again.
            obj.sunPosition = [-1,-1]; % we initialize the position of sun, to be used throughout the code
            obj.deg2radf = pi/180;
            obj.sunPattern = sunPattern;
            % create a video file reader
            try_ct = 0;
            obj.live = live;
            obj.vid_r = vid_r;
            obj.vid_n = vid_n;
            obj.vidObj = [];
            if(vid_r)
                obj.vidObj = VideoWriter(obj.vid_n);
                obj.vidObj.Quality = 100;
                obj.vidObj.FrameRate = 0.5;
                open(obj.vidObj);
            end
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
            obj.predict = predict; % whether we are in the tracking mode or display mode. can be 0 or 1
            obj.timeSeriesGHI = timeSeriesGHI; % the timeseries object corresponding to the GHI measurements
            obj.timeSeriesGHI_cs = timeSeriesGHI_cs;
            obj.timeSeriesPV = timeSeriesPV; % the timeseries object corresponding to the PV measurements
            obj.timeSeriesTemp = timeSeriesTemp; % the timeseries object corresponding to the temperature measurements
            obj.timeSeriesOcc = timeSeriesOcc; % the timeseries object corresponding to the occlusion measurements
            obj.timeSeriesCBH = timeSeriesCBH; % the timeseries object corresponding to the CBH measurements
            obj.irradianceValue = 0; % the variable that keeps the irradiance value corresponding to the current frame
            imageMask = ones(obj.calibModel.height,obj.calibModel.width);
            [X_l,Y_l] = meshgrid(1:obj.calibModel.height,1:obj.calibModel.width);
            obj.X = X_l;
            obj.Y = Y_l;
            P = [X_l(:)';Y_l(:)'];
            p = cam2world(P,obj.calibModel);
            p = obj.ecalibModel.R*p;
            maskTemp = p(3,:)<-0.15;
            indicesTemp = sub2ind([obj.calibModel.height,obj.calibModel.width],P(1,:),P(2,:));
            imageMask(indicesTemp) = maskTemp;
            imageMask = logical(repmat(imageMask,[1,1,3]));
            obj.imageMask = imageMask;
            obj.frameRate = 4;
            obj.timeHorizon = 25;
            obj.model3D = model3D;

            obj.occlusionPoints = [1:obj.timeHorizon];
            obj.occlusionTimes = [1:obj.timeHorizon]*obj.frameRate*3/(60);
            obj.occlusionTimeNums = obj.occlusionTimes/(60*24); %[1:obj.timeHorizon]*obj.frameRate*3/(60*60*24);
            obj.occlusionValues = zeros(1,obj.timeHorizon);
            obj.occlusionQualities = zeros(1,obj.timeHorizon);
            obj.attCs = 1;
            obj.attOcc = 2;
            [obj.shadowGrx,obj.shadowGry] = meshgrid([-5:0.5:5],[-5:0.5:5]);
            %[shadowGrx,shadowGry] = meshgrid([-2.5,-2,-1.5,0,1.5,2,2.5],[-2.5,-2,-1.5,0,1.5,2,2.5]);
            if(~obj.predict) % if the tracking mode is active initialize the needed display objects
                % create three video players, one to display the video,
                % one to display the foreground mask and objects, and one
                % to display the model.
                if(verbose)
                    obj.videoPlayer = VidPlayer([20, 500, 700, 400]); % for the real image
                    obj.maskPlayer = VidPlayer([740, 500, 700, 400]); % for the masked cloud-decision image
                    %obj.modPlayer = ModPlayer([740, 100, 700, 400],obj.ModPlayerOption); % for the model in real dimensions
                    obj.modPlayer = [];
                    obj.predPlayer = VidPlayer([740, 100, 700, 400]); % for the model in real dimensions
                    obj.mapPlayer = ModPlayer([20, 100, 700, 400],obj.ModPlayerOption,modelImage,obj.model3D.all); % for the model in real dimensions
                end
            else % othervise initialize the other objects
                % create two video players, one to display the video,
                % and one to display the irradiance values
                if(verbose)
                    obj.videoPlayer = VidPlayer([10, 550, 500, 400]); % for the real image
                    obj.maskPlayer = VidPlayer([550, 550, 500, 400]); % for the masked cloud-decision image
                    obj.predPlayer = VidPlayer([1090, 550, 500, 400]); % for the model in real dimensions
                    obj.plotPlayerTemp = ModPlayer([1090, 50, 500, 400],'plot'); % for the temperature-value plot
                    obj.plotPlayerPV = ModPlayer([550, 50, 500, 400],'plot'); % for the PV-value plot
                    obj.plotPlayerGHI = ModPlayer([10, 50, 500, 400],'plot'); % for the GHI-value plot
                    %obj.mapPlayer = ModPlayer([1750, 50, 1200, 900],obj.ModPlayerOption,modelImage,model3D.all);
                end
            end
        end
        
        function obj = readFrame(obj)
            % in this function we read the frame on the order from the folder
            
            % if the destination buffer do not exist, than we immediately
            % terminate
            if(exist(obj.reader) ~=7)
                disp('Input stream should have been closed... Terminating...');
                obj.frame = []; % indicator of an unsuccessful read trial
                return;
            end
            
            read_try = 1; % try to read a file for the first time
            if(obj.live) % if this is the live mode, meaning, we are going at the same paste as the camera
                obj.vid_ims = ls(obj.reader); % get the list of files in the folder
                obj.im_list = cellstr(obj.vid_ims(3:end,:)); % get the name list of those list of files
                while(length(obj.im_list{1}) == 0) % wait for an image to be written to the buffer
                    if(exist(obj.reader) ~=7) % if the buffer folder is deleted during the wait, this means the camera has been shut down, so terminate
                        disp('Input stream should have been closed... Terminating...');
                        obj.frame = []; % indicator of an unsuccessful read trial
                        return;
                    end
                    obj.vid_ims = ls(obj.reader); % get the list of files in the folder
                    obj.im_list = cellstr(obj.vid_ims(3:end,:)); % get the name list of those list of files
                    read_try = read_try+1; % increment the number of trials and go for another try
                    disp(char(['Trying to read for: ' num2str(read_try)])); % report the action to the user.
                    if(read_try > 8) % if a read is tried for more than 8 times, then terminate with an unsuccessfull read
                        disp('Input stream not available... Terminating...');
                        obj.frame = []; % unsuccessful read indicator
                        return;
                    end
                    pause(1); % wait 1 msec between two read trials, we can increase this time unit, but not so necessary
                end
                % wait until the name of the file in the buffer folder is
                % different from the previous read to make sure that the
                % observation frame is updated
                while(strcmp(char(strcat(obj.reader, '/', obj.im_list(1))),obj.frameName)==1)
                    pause(1);
                    obj.vid_ims = ls(obj.reader);
                    obj.im_list = cellstr(obj.vid_ims(3:end,:));
                end
                % if the code is here, this means that the new file is ready
                % to be read
                fid = fopen(char(strcat(obj.reader, '/', obj.im_list(1))), 'r'); % we try to read in case there is an operating system related lock on the file. This may happen if the camera is writing the
                % the file as we try to read it. if this application manages to open the file, than the file is locked by it and can be read.
                while (fid == -1) % if the camera is writing the file then apply the same trial based procedure until the write is over.
                    if(read_try > 8)
                        disp('Input stream not available... Terminating...');
                        obj.frame = [];
                        return;
                    end
                    disp(['Retrying from this one...   ' char(strcat(obj.reader, '/', obj.im_list(1)))]);
                    obj.vid_ims = ls(obj.reader);
                    obj.im_list = cellstr(obj.vid_ims(3:end,:));
                    fid = fopen(char(strcat(obj.reader, '/', obj.im_list(1))), 'r'); % try to lock the file once again
                    read_try = read_try + 1;
                    pause(0.5);
                end
                
                obj.frame = imread(char(strcat(obj.reader, '/', obj.im_list(1)))); % since the file is locked, it is safe to read it.
                fclose(fid); % since the file is read, the lock can be released
                
                obj.frameName = char(strcat(obj.reader, '/', obj.im_list(1))); % record the name of the file for the check operation on the up-to-date information of the file
                obj.frameCount = obj.frameCount + 1; % increase the frame count
            else
                % if it is not a live operation, directly read the file in
                % order
                if(length(obj.im_list) >= obj.frameCount) % if there are still some images to be read
                    obj.frameName = char(strcat(obj.reader, '/', obj.im_list(obj.frameCount)));
                    obj.frame = imread(obj.frameName);  % read the image directly, no need for the lock operations
                    obj.frameNameStr = strtok(char(obj.im_list(obj.frameCount)),'.'); % seperate the file type extension, get the date info
                    obj.frameDateVec = strread(obj.frameNameStr,'%d','delimiter','_'); % get the date vector elements, in integer, corresponding to day, month etc.
                    %obj.frameDateVec = frameDateVec;
                    
                    if(mod(obj.frameCount-obj.frameRate,50*obj.frameRate)<=obj.frameRate)
                        frameDateVecStr = strread(obj.frameNameStr,'%s','delimiter','_'); % get the date vector elements, in string, corresponding to day, month etc.
                        obj.frameRecordImageName = char(strcat(frameDateVecStr(1), '_',frameDateVecStr(2), '_',frameDateVecStr(3), ...
                            '_',frameDateVecStr(4), '_',frameDateVecStr(5))); % the name of the image or the gif file with the directory
                    end
                    obj.frameCount = obj.frameCount + obj.frameRate;
                else
                    obj.frame = []; % if no more frames in the recording, directly terminate
                    display('Done reading strored images... Terminating...');
                end
            end
            
            obj.frame(~obj.imageMask) = 0;  % we eliminate the parts that are angularly behind the camera observation point
        end
        
        function obj = prepareData(obj)
                        frameDateVec_l = obj.frameDateVec(1:6);
            frameDateVecStr = strread(obj.frameNameStr,'%s','delimiter','_'); % get the date vector elements, in string, corresponding to day, month etc.
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
            obj.frameDateStr = datestr(frameDateVec_l');
            queryDates = datestr([0;obj.occlusionTimeNums']+datenum(obj.frameDateStr));
            
            % read the irradiance value between the current frame date and minuteLag mins before
            %irradianceDB = obj.timeSeriesGHI.getsampleusingtime(frameDateStr,frameDateStr2); % query the values
            irradianceDB = resample(obj.timeSeriesGHI,queryDates);
            obj.irradianceValues = irradianceDB.Data(2:end);
            obj.irradianceTimes = obj.occlusionTimes';
            irradianceDB_cs = resample(obj.timeSeriesGHI_cs,queryDates);
            obj.irradianceValues_cs = irradianceDB_cs.Data(2:end);
            obj.irradianceTimes_cs = obj.occlusionTimes';
            if(isempty(obj.irradianceValues)) % if not found the associated entries, set to zero
                obj.irradianceValues = 0;
                obj.irradianceTimes = 0;
                obj.irradianceValues_cs = 0;
                obj.irradianceTimes_cs = 0;
            end
            
            % read the power value between the current frame date and minuteLag mins before.
            powerDB = resample(obj.timeSeriesPV,queryDates);
            obj.powerValues = powerDB.Data(2:end); % get the values into buffer
            obj.powerTimes = obj.occlusionTimes';
            if(isempty(obj.powerValues)) % if not found the associated entries, set to zero
                obj.powerValues = 0;
                obj.powerTimes = 0;
            end
            
            % read the temperature value between the current frame date and minuteLag mins before.
            tempDB = resample(obj.timeSeriesTemp,queryDates);
            obj.tempValues = tempDB.Data(2:end); % get the values into buffer
            obj.tempTimes = obj.occlusionTimes';
            if(isempty(obj.tempValues)) % if not found the associated entries, set to zero
                obj.tempValues = 0;
                obj.tempTimes = 0;
            end
            
            cbhDB = resample(obj.timeSeriesCBH, obj.frameDateStr); % query the values
            obj.cbh =  cbhDB.Data(1);
            
            occDB = resample(obj.timeSeriesOcc,queryDates);
            obj.occValues = occDB.Data(2:end,1)>10 | occDB.Data(2:end,2)<-10;
            obj.occTimes = obj.occlusionTimes';
            if(isempty(obj.occValues)) % if not found the associated entries, set to zero
                obj.occValues = 0;
                obj.occTimes = 0;
            end
            
            obj.occState = occDB.Data(1,1)>10 | occDB.Data(1,2)<-10;
            obj.clearnessState = ~(occDB.Data(1,1)>10 | occDB.Data(1,2)<-10) & (occDB.Data(1,1)<-10 | occDB.Data(1,2)>=0);

            obj.irrCs = irradianceDB_cs.Data(1);
            obj.irrMs = irradianceDB.Data(1);
            
            if(obj.occState)
                if(obj.attOcc == 2)
                    obj.attOcc = obj.irrCs/obj.irrMs;
                else
                    obj.attOcc = 0.7*obj.attOcc + 0.3*obj.irrCs/obj.irrMs;
                    if(obj.attOcc < obj.attCs)
                        obj.attOcc = 1.5;
                    end
                    %obj.attOcc = obj.irradianceValues_cs(1)/obj.irradianceValues(1);
                end
            elseif(obj.clearnessState)
                if(obj.attCs == 1)
                    obj.attCs = obj.irrCs/obj.irrMs;
                else
                    obj.attCs = 0.9*obj.attCs + 0.1*obj.irrCs/obj.irrMs;
                    if(obj.attCs > 2)
                        obj.attCs = 1.2;
                        %elseif(obj.attCs < 1)
                        %    obj.attCs = 1.01;
                    end
                    %obj.attCs = obj.irradianceValues_cs(1)/obj.irradianceValues(1);
                end
            end
%            if(obj.occState)
%                 if(obj.attOcc == 2)
%                     obj.attOcc = obj.irrCs/obj.irrMs;
%                 else
%                     obj.attOcc = 0.7*obj.attOcc + 0.3*obj.irrCs/obj.irrMs;
%                     if(obj.attOcc < 1)
%                         obj.attOcc = 1.7;
%                     end
%                     obj.attOcc = obj.irradianceValues_cs(1)/obj.irradianceValues(1);
%                 end
%             else
%                 if(obj.attCs == 1)
%                     obj.attCs = obj.irrCs/obj.irrMs;
%                 else
%                     obj.attCs = 0.7*obj.attCs + 0.3*obj.irrCs/obj.irrMs;
%                     if(obj.attCs > 2)
%                         obj.attCs = 1.3;
%                     end
%                     obj.attCs = obj.irradianceValues_cs(1)/obj.irradianceValues(1);
%                 end
%             end
            
            obj.irradianceValues_est = obj.irradianceValues_cs;
        end
        
        function obj = detectObjects(obj)
           
            %contours = cv.findContours(mask,'Mode','List'); % contour extraction using openCV method
            obj.contours = cv.findContours(obj.mask,'Mode','External'); % contour extraction using openCV method
            % areas of objects are approximately computed here using several methods
            %areas = cell2mat(cellfun(@(x) size(cell2mat(x(:)),1),contours,'UniformOutput',false))';
            %areas = cell2mat(cellfun(@(x) prod(max(cell2mat(x(:)),[],1)-min(cell2mat(x(:)),[],1)) ,contours,'UniformOutput',false))';
            obj.areas = cell2mat(cellfun(@(x) cv.contourArea(x(:)) ,obj.contours,'UniformOutput',false))';
            
            % elminate the objects with contour smaller than 350.
            consider = (obj.areas > 300);% & (areas < 100500); % if the area of an object is less than 350 pixels, we discard that
            
            % take only the objects that we want to consider
            obj.contours = obj.contours(consider);
            obj.areas = obj.areas(consider);
            %+ones(size(cell2mat(x(:))))
            obj.contoursProj = cellfun(@(x) num2cell(cloudTrackingCV.getImage2Real(cell2mat(x(:))+1,obj.calibModel,obj.ecalibModel,obj.cbh),2)',obj.contours,'UniformOutput',false);
            obj.areasProj = cell2mat(cellfun(@(x) cv.contourArea(x(:)) ,obj.contoursProj,'UniformOutput',false))';
            
            % different methods for determining the point that represent the
            % cloud objets first method uses the mean of the contours the
            % second one uses the mid point of the bounding rectangle
            %obj.centroids = cellfun(@(x) mean(cell2mat(x(:)),1)+1,obj.contours,'UniformOutput',false);
            %obj.centroidsProj = cellfun(@(x) mean(cell2mat(x(:)),1),obj.contoursProj,'UniformOutput',false);
            obj.centroids = cellfun(@(x) (min(cell2mat(x(:))+1,[],1)+max(cell2mat(x(:))+1,[],1))/2,obj.contours,'UniformOutput',false);
            obj.centroidsProj = cellfun(@(x) (min(cell2mat(x(:)),[],1)+max(cell2mat(x(:)),[],1))/2,obj.contoursProj,'UniformOutput',false);
            
            obj.bboxes = cellfun(@(x) [min(cell2mat(x(:))+1,[],1),max(cell2mat(x(:))+1,[],1)-min(cell2mat(x(:))+1,[],1)] ,obj.contours,'UniformOutput',false);
            obj.bboxesProj = cellfun(@(x) [min(cell2mat(x(:)),[],1),max(cell2mat(x(:)),[],1)-min(cell2mat(x(:)),[],1)] ,obj.contoursProj,'UniformOutput',false);
            
            obj.centroids = vertcat(obj.centroids{:});
            obj.centroidsProj = vertcat(obj.centroidsProj{:});
            obj.bboxes = vertcat(obj.bboxes{:}); % put all the boxes in a 2D array, so that the annotation algorithm can display them
            obj.bboxesProj = vertcat(obj.bboxesProj{:});
            
        end
        
        function obj = getBinaryImage(obj, type)
            r_b = double(obj.frame(:,:,1))./double(obj.frame(:,:,3)); % R/B thresholding applied
            rgb = double(obj.frame(:,:,1))+double(obj.frame(:,:,2))+double(obj.frame(:,:,3)); % to distinguish other objects like buildings from clouds, a test on brightness
            obj.mask = r_b > 0.72 & rgb > 100;             %% brightness auto thing is needed because of this ..
            
            sunP = obj.sunPatch;
            g_b = double(obj.frame(:,:,2))./double(obj.frame(:,:,3));
            mask2 = g_b > 0.8;
            %mask2 = (r_b - g_b)>-0.01;
            obj.mask(sunP) = mask2(sunP);
            
            if(strcmp(type,'edge'))
                H = [1,1,1;1,-8,1;1,1,1];
                obj.mask = imfilter(obj.mask,H,'replicate');
            end
            obj.mask = logical(cv.medianBlur(obj.mask, 'KSize', 13));
        end
        
        function obj = getSunPositionReal(obj)
            DN = datenum(obj.frameDateVec(1:6)'); % date number corresponding to the current frame
            Time = pvl_maketimestruct(DN, obj.model3D.UTC); % time structure corresponding to the current frame
            [SunAz, ~, ApparentSunEl, ~]=pvl_ephemeris(Time, obj.model3D.Location); % sun angles corresponding to the current frame time and the camera location
            zenith = 90-ApparentSunEl; % sun zenith angle
            azimuth = SunAz; % sun azimuth angle
            [x,y,z] = sph2cart((90-azimuth)*obj.deg2radf,(90-zenith)*obj.deg2radf,1); % direction of the sun converted from the angles to spherical coordinates
            sunPositionR = [x;y;-z];
            sunPositionReal = world2cam(obj.ecalibModel.Rinv*sunPositionR,obj.calibModel); % direction of the sun converted from global to image pixel coordinates
            obj.sunPosition = [sunPositionReal(2);sunPositionReal(1)];
            [x_sd,y_sd,~] = sph2cart((90-azimuth)*obj.deg2radf,(90-zenith)*obj.deg2radf,obj.cbh/cos(zenith*obj.deg2radf)); % global position of the sun is computed assuming a certain height
            obj.shadowDisp = [-x_sd,-y_sd]; % the shadow displacement depends on the global position of the sun and the assumed cbh 
        end
    
        function obj = getSunPosition(obj)
            res = cv.matchTemplate(obj.frame,obj.sunPattern);
            if(min(min(res))>1e8)
                min(min(res))
            end
            [y,x]=ind2sub(size(res),find(res==min(min(res))));
            
            s_s = size(obj.sunPattern);
            x_s = s_s(2);
            y_s = s_s(1);
            if(size(x)>(size(res,1)/3))
                obj.sunPosition = [-1,-1];
                return;
            end
            obj.sunPosition = [x,y]+([x_s,y_s]./2);
        end
        
        function obj = getSunPatch(obj)
            obj.sunPatch = sqrt((obj.X-obj.sunPosition(1)).^2+(obj.Y-obj.sunPosition(2)).^2)<=200;
        end
        
        function obj = eliminateSun(obj)
            x = obj.sunPosition(1);
            y = obj.sunPosition(2);
            obj.frame = cv.circle(obj.frame,[x,y],80,'Color',[0,0,0],'Thickness',-1);
        end
        
        function obj = getShadowDisplacement(obj)
            DN = datenum(obj.frameDateVec(1:6)');
            Time = pvl_maketimestruct(DN, obj.model3D.UTC);
            [SunAz, ~, ApparentSunEl, ~]=pvl_ephemeris(Time, obj.model3D.Location);
            zenith = 90-ApparentSunEl;
            azimuth = SunAz;
            [x,y,~] = sph2cart((90-azimuth)*obj.deg2radf,(90-zenith)*obj.deg2radf,obj.cbh/cos(zenith*obj.deg2radf));
            obj.shadowDisp = [-x,-y];
        end
        
        function obj = getShadowDisplacement_DN(obj,DN)
            Time = pvl_maketimestruct(DN, obj.model3D.UTC);
            [SunAz, ~, ApparentSunEl, ~]=pvl_ephemeris(Time, obj.model3D.Location);
            zenith = 90-ApparentSunEl;
            azimuth = SunAz;
            [x,y,~] = sph2cart((90-azimuth)*obj.deg2radf,(90-zenith)*obj.deg2radf,obj.cbh/cos(zenith*obj.deg2radf));
            obj.shadowDisp_d = [-x,-y];
        end
        
        function obj = predictNewLocationsOfTracksIntoFuture(obj)
            %tracksPred = [tracks,tracksLost];
            obj.tracksPred = obj.tracks;
            obj.predCloudContours = [];
            obj.predShadowContours = [];
            %occlusion = zeros(1,timeSteps);
            obj.occlusionValues = [obj.occlusionValues(2:end),0];
            obj.occlusionQualities = [obj.occlusionQualities(2:end),0];
            for predNum = 1:obj.timeHorizon
                occ = 0;
                upInd = find(obj.occlusionPoints==predNum,1);
                for i = 1:length(obj.tracksPred)
                    % to be able to used them for image annotation and sanity
                    % checking, we provide the bounding boxes in an array format
                    % where one bounding box is on one row of the matrix
                    bbox = obj.tracksPred(i).bbox;
                    bboxProj = obj.tracksPred(i).bboxProj;
                    
                    % predict the current location of the track
                    [obj.tracksPred(i).kalmanFilter,predictedCentroid,~] = obj.tracksPred(i).kalmanFilter.Estimate();
                    obj.tracksPred(i).centroidFilter = cloudTrackingCV.filterAdd(predictedCentroid, obj.tracksPred(i).centroidFilter);  % the extra filters are not operating right now, they are put in case needed.
                    obj.tracksPred(i).predCentroid = obj.tracksPred(i).centroidFilter.value;
                    
                    % projected values from the kalman filter.
                    [obj.tracksPred(i).kalmanFilterProj,predictedCentroidProj,~] = obj.tracksPred(i).kalmanFilterProj.Estimate();
                    obj.tracksPred(i).centroidProjFilter = cloudTrackingCV.filterAdd(predictedCentroidProj, obj.tracksPred(i).centroidProjFilter);  % the extra filters are not operating right now, they are put in case needed.
                    obj.tracksPred(i).predCentroidProj = obj.tracksPred(i).centroidProjFilter.value;
                    
                    % here we update the state variables of the current track
                    % according to different types of the kalman filter
                    switch obj.trackOption
                        case {'force', 'acceleration'}
                            obj.tracksPred(i).predSpeed = [obj.tracksPred(i).kalmanFilter.State(2),obj.tracksPred(i).kalmanFilter.State(5)];
                        case 'speed'
                            obj.tracksPred(i).predSpeed = [obj.tracksPred(i).kalmanFilter.State(2),obj.tracksPred(i).kalmanFilter.State(4)];
                        otherwise
                    end
                    
                    % here to update bounding boxes, we get a local copy on the
                    % location of tracks
                    predictedCentroid = obj.tracksPred(i).predCentroid;
                    
                    % shift the bounding box so that its center is at
                    % the predicted location
                    predictedCentroid = int32(predictedCentroid(1:2)) - (int32(bbox(3:4)) / 2);
                    
                    % reform  the bounding boxes with their new locations
                    obj.tracksPred(i).bbox = double([predictedCentroid(1:2), bbox(3:4)]);
                    
                    % here to update bounding boxes (on the projected plane), we get a local copy on the
                    % location of tracks
                    predictedCentroidProj = obj.tracksPred(i).predCentroidProj;
                    
                    % reform  the bounding boxes with their new locations
                    obj.tracksPred(i).bboxProj = double([predictedCentroidProj(1:2), bboxProj(3:4)]);
                    
                    if(~isempty(upInd))
                        obj = getShadowDisplacement_DN(obj,datenum(obj.frameDateVec(1:6)')+ predNum*3*obj.frameRate/(60*60*24));
                        cloudContours = [obj.tracksPred(i).contourProj+repmat(obj.tracksPred(i).predCentroidProj,size(obj.tracksPred(i).contourProj,1),1),obj.cbh*ones(size(obj.tracksPred(i).contourProj,1),1)];
                        shadowContours = [obj.tracksPred(i).contourProj+repmat(obj.tracksPred(i).predCentroidProj+obj.shadowDisp_d,size(obj.tracksPred(i).contourProj,1),1),zeros(size(obj.tracksPred(i).contourProj,1),1)];
                        obj.predCloudContours{1,i} = num2cell(cloudContours,2)';
                        obj.predShadowContours{1,i} = num2cell(shadowContours,2)';
                        [in,on] = inpolygon(obj.model3D.Sen_xs+obj.shadowGrx,obj.model3D.Sen_ys+obj.shadowGry,shadowContours(:,1),shadowContours(:,2));
                        %[in,on] = inpolygon(model3D.Sen_xs,model3D.Sen_ys,shadowContours(:,1),shadowContours(:,2));
                        occRes = in | on;
                        occ = occ + sum(occRes(:));
                        %occlusion(upInd) = occlusion(upInd) | sum(occRes(:))>0;
                    end
                    %occlusion(predNum) =  double(sum(occRes(:)))/length(occRes(:));
                end
                if(~isempty(upInd) && occ>0)
                    if(obj.occlusionValues(predNum)==0)
                        obj.occlusionQualities(predNum) = occ;
                    end
                    obj.occlusionValues(predNum) = obj.occlusionValues(predNum) | occ>0;
                else
                    if(obj.occlusionQualities(predNum) < length(obj.shadowGry(:)))
                        obj.occlusionValues(predNum) = 0.96*obj.occlusionValues(predNum) + 0.04*occ;
                    end
                end
            end
            obj.irradianceValues_est = obj.irradianceValues_est ./ (obj.occlusionValues'*obj.attOcc + (1-obj.occlusionValues')*obj.attCs);
        end
        
        function obj = getPredictedImage(obj)
            obj.predFrame = zeros(size(obj.frame));
            predContours=[];
            iterator = 1;
            for i = 1:length(obj.tracksPred)
                if(obj.tracksPred(i).consecutiveInvisibleCount<1)
                    predContours{1,iterator} = num2cell(cloudTrackingCV.getReal2Image(obj.tracksPred(i).contourProj+repmat(obj.tracksPred(i).predCentroidProj,size(obj.tracksPred(i).contourProj,1),1),obj.calibModel,obj.ecalibModel,obj.cbh)-1,2)';
                    iterator = iterator+1;
                end
            end
            if(length(obj.tracksPred)~=0)
                obj.predFrame = cv.drawContours(obj.predFrame,predContours,'Thickness',-1,'LineType',4);
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
        
        function obj = predictNewLocationsOfTracks(obj)
            for i = 1:length(obj.tracks)
                obj.tracks(i).cbh = obj.cbh;
                % to be able to used them for image annotation and sanity
                % checking, we provide the bounding boxes in an array format
                % where one bounding box is on one row of the matrix
                bbox = obj.tracks(i).bbox;
                bboxProj = obj.tracks(i).bboxProj;
                
                % predict the current location of the track
                [obj.tracks(i).kalmanFilter,predictedCentroid,~] = obj.tracks(i).kalmanFilter.Estimate();
                obj.tracks(i).centroidFilter = cloudTrackingCV.filterAdd(predictedCentroid, obj.tracks(i).centroidFilter);  % the extra filters are not operating right now, they are put in case needed.
                obj.tracks(i).predCentroid = obj.tracks(i).centroidFilter.value;
                
                % projected values from the kalman filter.
                [obj.tracks(i).kalmanFilterProj,predictedCentroidProj,~] = obj.tracks(i).kalmanFilterProj.Estimate();
                obj.tracks(i).centroidProjFilter = cloudTrackingCV.filterAdd(predictedCentroidProj, obj.tracks(i).centroidProjFilter);  % the extra filters are not operating right now, they are put in case needed.
                obj.tracks(i).predCentroidProj = obj.tracks(i).centroidProjFilter.value;
                
                %tracks(i).predCentroid = getReal2Image(tracks(i).predCentroidProj,obj.calibModel,obj.ecalibModel,cbh);
                % here we update the state variables of the current track
                % according to different types of the kalman filter
                switch obj.trackOption
                    case {'force', 'acceleration'}
                        obj.tracks(i).predSpeed = [obj.tracks(i).kalmanFilter.State(2),obj.tracks(i).kalmanFilter.State(5)];
                    case 'speed'
                        obj.tracks(i).predSpeed = [obj.tracks(i).kalmanFilter.State(2),obj.tracks(i).kalmanFilter.State(4)];
                    otherwise
                end
                
                % here to update bounding boxes, we get a local copy on the
                % location of tracks
                predictedCentroid = obj.tracks(i).predCentroid;
                
                % shift the bounding box so that its center is at
                % the predicted location
                predictedCentroid = int32(predictedCentroid(1:2)) - (int32(bbox(3:4)) / 2);
                
                % reform  the bounding boxes with their new locations
                obj.tracks(i).bbox = double([predictedCentroid(1:2), bbox(3:4)]);
                
                % here to update bounding boxes (on the projected plane), we get a local copy on the
                % location of tracks
                predictedCentroidProj = obj.tracks(i).predCentroidProj;
                
                % shift the bounding box so that its center is at
                % the predicted location
                predictedCentroidProj = predictedCentroidProj(1:2) - bboxProj(3:4) / 2.0;
                
                % reform  the bounding boxes with their new locations
                obj.tracks(i).bboxProj = double([predictedCentroidProj(1:2), bboxProj(3:4)]);
            end
        end
        
        function obj = detectionToTrackAssignment(obj)
            % in this function we find the assignment between the existing
            % tracks and the observations coming from a certain frame. We do
            % the assignment on the image plane, so the numbers correspond to
            % pixels mostly.
            
            nTracks = length(obj.tracks);
            nDetections = size(obj.centroidsProj, 1);
            % we initialize the assignment matrix and all other related
            % variables first
            if(nTracks == 0 || size(obj.centroidsProj,1) == 0)
                obj.assignments = [];
                obj.unassignedTracks = [];
                obj.unassignedDetections = [1:nDetections];
                return;
            end
            
            % compute the cost of assigning each detection to each track
            % here we want to obtain (x_i - y_j)^2 where x_i is the position of
            % a track and y_j is the position of an observation. for this
            % purpose we first get the x_i^2, then we get the y_j^2, after we
            % get the x_i*y_j. after those it all comes down to calculating
            % x_i^2 - 2*x_i*y_j + y_j^2 using the matrices.
            cent_tracks = cat(1,obj.tracks(:).predCentroidProj);
            match_mat1 = dot(cent_tracks,cent_tracks,2); % square of tracks x_i^2
            match_mat2 = dot(obj.centroidsProj,obj.centroidsProj,2); % square of detections y_j^2
            match_mat11 = repmat(match_mat1,1,length(match_mat2)); % get x_i^2 into matrix
            match_mat22 = repmat(match_mat2',length(match_mat1),1); % get y_j^2 into matrix
            match_mat12 = cent_tracks*obj.centroidsProj'; % get x_i*y_j into matrix
            
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
            obj.assignments = [id_tracks(matching~=0),matching(matching~=0)]; % here we get the assigned tracks
            obj.unassignedTracks = setdiff(id_tracks,obj.assignments(:,1)); % get the tracks that are unassigned
            obj.unassignedDetections = setdiff([1:nDetections],obj.assignments(:,2)); % get the observations that are not assigned
        end
        
        function obj = detectionToTrackAssignmentDetailed(obj)
            % in this function we find the assignment between the existing
            % tracks and the observations coming from a certain frame. We do
            % the assignment on the image plane, so the numbers correspond to
            % pixels mostly.
            
            nTracks = length(obj.tracks);
            nDetections = size(obj.centroidsProj, 1);
            % we initialize the assignment matrix and all other related
            % variables first
            if(nTracks == 0 || size(obj.centroidsProj,1) == 0)
                obj.assignments = [];
                obj.unassignedTracks = [];
                obj.unassignedDetections = [1:nDetections];
                return;
            end
            
            % compute the cost of assigning each detection to each track
            % here we want to obtain (x_i - y_j)^2 where x_i is the position of
            % a track and y_j is the position of an observation. for this
            % purpose we first get the x_i^2, then we get the y_j^2, after we
            % get the x_i*y_j. after those it all comes down to calculating
            % x_i^2 - 2*x_i*y_j + y_j^2 using the matrices.
            cent_tracks = [cat(1,obj.tracks(:).predCentroidProj),cat(1,obj.tracks(:).bboxProj)];
            cent_tracks(:,5:6) = cent_tracks(:,3:4)+cent_tracks(:,5:6);
            cent_proj = [obj.centroidsProj,obj.bboxesProj];
            cent_proj(:,5:6) = cent_proj(:,3:4)+cent_proj(:,5:6);
            normalization_factor = max([cent_tracks;cent_proj],[],1) - min([cent_tracks;cent_proj],[],1);
            cent_tracks_nrm = cent_tracks*diag(normalization_factor);
            cent_proj_nrm = cent_proj*diag(normalization_factor);
            match_mat1 = dot(cent_tracks_nrm,cent_tracks_nrm,2); % square of tracks x_i^2
            match_mat2 = dot(cent_proj_nrm,cent_proj_nrm,2); % square of detections y_j^2
            match_mat11 = repmat(match_mat1,1,length(match_mat2)); % get x_i^2 into matrix
            match_mat22 = repmat(match_mat2',length(match_mat1),1); % get y_j^2 into matrix
            match_mat12 = cent_tracks_nrm*cent_proj_nrm'; % get x_i*y_j into matrix
            
            % we construct the cost matrix here. entry i j refers to the cost
            % of assigning the i'th track to the j'th observation
            %cost = match_mat11 - 2*match_mat12 + match_mat22; % get x_i^2 - 2*x_i*y_j + y_j^2 using those matrices
            cost = 1 - match_mat12;
            %cost_th = 9*(obj.cbh^2/300^2)*1e4;
            %cost(cost>cost_th) = inf; % here we do a sanity check on the costs. if the cost of assigning is greater than 2500 meaning
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
            obj.assignments = [id_tracks(matching~=0),matching(matching~=0)]; % here we get the assigned tracks
            obj.unassignedTracks = setdiff(id_tracks,obj.assignments(:,1)); % get the tracks that are unassigned
            obj.unassignedDetections = setdiff([1:nDetections],obj.assignments(:,2)); % get the observations that are not assigned
        end
        
        function obj = detectionToTrackAssignmentIntersect(obj)
            % in this function we find the assignment between the existing
            % tracks and the observations coming from a certain frame. We do
            % the assignment on the image plane, so the numbers correspond to
            % pixels mostly.
            
            nTracks = length(obj.tracks);
            nDetections = size(obj.centroidsProj, 1);
            % we initialize the assignment matrix and all other related
            % variables first
            if(nTracks == 0 || size(obj.centroidsProj,1) == 0)
                obj.assignments = [];
                obj.unassignedTracks = [];
                obj.unassignedDetectionsToTracks = [];
                obj.unassignedTracksToDetections = [];
                obj.unassignedDetections = [1:nDetections];
                return;
            end
            
            % compute the cost of assigning each detection to each track
            % here we want to obtain (x_i - y_j)^2 where x_i is the position of
            % a track and y_j is the position of an observation. for this
            % purpose we first get the x_i^2, then we get the y_j^2, after we
            % get the x_i*y_j. after those it all comes down to calculating
            % x_i^2 - 2*x_i*y_j + y_j^2 using the matrices.
            cent_tracks = cat(1,obj.tracks(:).predCentroidProj);
            match_mat1 = dot(cent_tracks,cent_tracks,2); % square of tracks x_i^2
            match_mat2 = dot(obj.centroidsProj,obj.centroidsProj,2); % square of detections y_j^2
            match_mat11 = repmat(match_mat1,1,length(match_mat2)); % get x_i^2 into matrix
            match_mat22 = repmat(match_mat2',length(match_mat1),1); % get y_j^2 into matrix
            match_mat12 = cent_tracks*obj.centroidsProj'; % get x_i*y_j into matrix
            
            bbox_tracks = cat(1,obj.tracks(:).bboxProj);
            bbox_tracks(:,3:4) = bbox_tracks(:,1:2)+bbox_tracks(:,3:4);
            bbox_proj = obj.bboxesProj;
            bbox_proj(:,3:4) = bbox_proj(:,1:2)+bbox_proj(:,3:4);
            
            int_mat11 = repmat(bbox_tracks(:,1),1,length(bbox_proj(:,1)));
            int_mat12 = repmat(bbox_tracks(:,2),1,length(bbox_proj(:,2)));
            int_mat13 = repmat(bbox_tracks(:,3),1,length(bbox_proj(:,3)));
            int_mat14 = repmat(bbox_tracks(:,4),1,length(bbox_proj(:,4)));
            
            int_mat21 = repmat(bbox_proj(:,1)',length(bbox_tracks(:,1)),1);
            int_mat22 = repmat(bbox_proj(:,2)',length(bbox_tracks(:,2)),1);
            int_mat23 = repmat(bbox_proj(:,3)',length(bbox_tracks(:,3)),1);
            int_mat24 = repmat(bbox_proj(:,4)',length(bbox_tracks(:,4)),1);
            
            int_mat = (int_mat11 < int_mat23) & (int_mat12 < int_mat24) & (int_mat13 > int_mat21) & (int_mat14 > int_mat22);
            
            area_tracks = cat(1,obj.tracks(:).areaProj);
            area_mat1 = repmat(area_tracks,1,length(obj.areasProj));
            area_mat2 = repmat(obj.areasProj',length(area_tracks),1);
            area_mat = area_mat2./area_mat1;
            area_mat = (area_mat>0.2) & (area_mat<5);
            % we construct the cost matrix here. entry i j refers to the cost
            % of assigning the i'th track to the j'th observation
            cost = match_mat11 - 2*match_mat12 + match_mat22; % get x_i^2 - 2*x_i*y_j + y_j^2 using those matrices
            %cost_th = (obj.cbh^2/300^2)*1e4;
            %cost(cost>cost_th) = inf; % here we do a sanity check on the costs. if the cost of assigning is greater than 2500 meaning
            cost(~(int_mat & area_mat)) = inf;
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
            obj.assignments = [id_tracks(matching~=0),matching(matching~=0)]; % here we get the assigned tracks
            if(~isempty(obj.assignments))
                obj.unassignedTracks = setdiff(id_tracks,obj.assignments(:,1)); % get the tracks that are unassigned
                obj.unassignedDetections = setdiff([1:nDetections],obj.assignments(:,2))'; % get the observations that are not assigned
            else
                obj.unassignedTracks = id_tracks; % get the tracks that are unassigned
                obj.unassignedDetections = [1:nDetections]'; % get the observations that are not assigned
            end
            [unDeCosts,unDeTracks] = min(cost(:,obj.unassignedDetections),[],1);
            if(~isempty(unDeCosts))
                obj.unassignedDetectionsToTracks = [unDeTracks(unDeCosts~=inf)',obj.unassignedDetections(unDeCosts~=inf)];
            else
                obj.unassignedDetectionsToTracks = [];
            end
            [unTCosts,unTDetections] = min(cost(obj.unassignedTracks,:),[],2);
            if(~isempty(unTCosts))
                obj.unassignedTracksToDetections = [obj.unassignedTracks(unTCosts~=inf),unTDetections(unTCosts~=inf)];
            else
                obj.unassignedTracksToDetections = [];
            end
        end

        function obj = detectionToTrackAssignmentGaussian(obj)
            % in this function we find the assignment between the existing
            % tracks and the observations coming from a certain frame. We do
            % the assignment on the image plane, so the numbers correspond to
            % pixels mostly.
            
            nTracks = length(obj.tracks);
            nDetections = size(obj.centroidsProj, 1);
            % we initialize the assignment matrix and all other related
            % variables first
            if(nTracks == 0 || size(obj.centroidsProj,1) == 0)
                obj.assignments = [];
                obj.unassignedTracks = [];
                obj.unassignedDetections = [1:nDetections];
                return;
            end
            
            
            for i = 1:nTracks
                P = obj.tracks(i).kalmanFilterProj.P;
                P_s = [P(1,1),P(1,3),0;P(3,1),P(3,3),0;0,0,1e6];
                cost(i,:) = 1-mvnpdf([obj.centroidsProj,obj.areas],[obj.tracks(i).predCentroidProj,obj.tracks(i).area],P_s)';
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
            obj.assignments = [id_tracks(matching~=0),matching(matching~=0)]; % here we get the assigned tracks
            obj.unassignedTracks = setdiff(id_tracks,obj.assignments(:,1)); % get the tracks that are unassigned
            obj.unassignedDetections = setdiff([1:nDetections],obj.assignments(:,2)); % get the observations that are not assigned
        end
        
        function obj = updateAssignedTracks(obj)
            % here we update then statistical parameters corresponding to
            % assgined tracks along with the kalman filters of those tracks
            % with the observartions associated
            numAssignedTracks = size(obj.assignments, 1);
            for i = 1:numAssignedTracks % for each assignment
                trackIdx = obj.assignments(i, 1); % get the track that is assigned
                detectionIdx = obj.assignments(i, 2); % get the detection id that is assigned
                centroid = obj.centroids(detectionIdx, :); % get the observation correspoding to the position of the detected object, on image plane
                centroidProj = obj.centroidsProj(detectionIdx, :); % get the observation correspoding to the position of the detected object, on projected plane
                bbox = obj.bboxes(detectionIdx, :); % get the observation correspoding to the bounding box of the detected object, on image plane
                bboxProj = obj.bboxesProj(detectionIdx, :); % get the observation correspoding to the bounding box of the detected object, on projected plane
                area = obj.areas(detectionIdx, :); % get the observation correspoding to the area of the detected object, on projected plane
                areaProj = obj.areasProj(detectionIdx, :); % get the observation correspoding to the area of the detected object, on projected plane
                contour = vertcat(obj.contours{1,detectionIdx}{1,:});
                contour = contour - repmat(centroid,size(contour,1),1);
                contourProj = vertcat(obj.contoursProj{1,detectionIdx}{1,:});
                contourProj = contourProj - repmat(centroidProj,size(contourProj,1),1);
                
                obj.tracks(trackIdx).centroidObs = centroid;
                obj.tracks(trackIdx).centroidProjObs = centroidProj;
                obj.tracks(trackIdx).contour = contour;
                obj.tracks(trackIdx).contourProj = contourProj;
                
                areaRatio = area/obj.tracks(trackIdx).area;
                if( areaRatio < 0.6 && areaRatio > 1.67 )
                	obj.tracks(trackIdx).kalmanFilter.State(obj.tracks(trackIdx).kalmanFilter.state_to_x) = centroid';
                    obj.tracks(trackIdx).kalmanFilterProj.State(obj.tracks(trackIdx).kalmanFilterProj.state_to_x) = centroidProj';
                end
                
                obj.tracks(trackIdx).prevCentroid = obj.tracks(trackIdx).centroid; % get a recording on the previous position of the detected object, on image plane
                obj.tracks(trackIdx).prevCentroidProj = obj.tracks(trackIdx).centroidProj; % get a recording on the previous position of the detected object, on projected plane
                
                % correct the estimate of the object's location
                % using the new detection
                
                % filter out the position observation
                obj.tracks(trackIdx).positionFilter = cloudTrackingCV.filterAdd(centroid,obj.tracks(trackIdx).positionFilter);
                pos_t = obj.tracks(trackIdx).positionFilter.value;
                
                % filter out the area observation
                obj.tracks(trackIdx).massFilter = cloudTrackingCV.filterAdd(area,obj.tracks(trackIdx).massFilter);
                obj.tracks(trackIdx).area = obj.tracks(trackIdx).massFilter.value;
                
                obj.tracks(trackIdx).massFilterProj = cloudTrackingCV.filterAdd(area,obj.tracks(trackIdx).massFilterProj);
                obj.tracks(trackIdx).areaProj = obj.tracks(trackIdx).massFilterProj.value;
                
                % fiter out the centroid information
                obj.tracks(trackIdx).centroidFilter = cloudTrackingCV.filterAdd(centroid,obj.tracks(trackIdx).centroidFilter);
                obj.tracks(trackIdx).centroid = obj.tracks(trackIdx).centroidFilter.value;
                
                % fiter out the projected centroid information
                obj.tracks(trackIdx).centroidProjFilter = cloudTrackingCV.filterAdd(centroidProj,obj.tracks(trackIdx).centroidProjFilter);
                obj.tracks(trackIdx).centroidProj = obj.tracks(trackIdx).centroidProjFilter.value;
                
                % get the observed speed wrt the previous step using a first
                % order differential approximation
                obj.tracks(trackIdx).obsSpeed = obj.tracks(trackIdx).centroid - obj.tracks(trackIdx).prevCentroid;
                
                % update the kalman filter corresponding to the centroid, on
                % the image plane
                obj.tracks(trackIdx).kalmanFilter = obj.tracks(trackIdx).kalmanFilter.Update(obj.tracks(trackIdx).centroid',double(obj.tracks(trackIdx).area));
                obj.tracks(trackIdx).centroid = obj.tracks(trackIdx).kalmanFilter.x';
                
                % update the kalman filter corresponding to the centroid, on
                % the projected
                obj.tracks(trackIdx).kalmanFilterProj = obj.tracks(trackIdx).kalmanFilterProj.Update(obj.tracks(trackIdx).centroidProj',double(obj.tracks(trackIdx).area));
                obj.tracks(trackIdx).centroidProj = obj.tracks(trackIdx).kalmanFilterProj.x';
                
                % getting the speed estimate from the kalman depending on the type of kalman used...
                switch obj.trackOption
                    case {'force','acceleration'}
                        speed = [obj.tracks(trackIdx).kalmanFilter.State(2), obj.tracks(trackIdx).kalmanFilter.State(5)];
                    case 'speed'
                        speed = [obj.tracks(trackIdx).kalmanFilter.State(2), obj.tracks(trackIdx).kalmanFilter.State(4)];
                    otherwise
                end
                
                % filter out the speed information
                obj.tracks(trackIdx).speedFilter = cloudTrackingCV.filterAdd(speed, obj.tracks(trackIdx).speedFilter);
                speed = obj.tracks(trackIdx).speedFilter.value;% * camParams.conversionFactor;
                
                % for some reason, in case needed, get the x and y speeds
                obj.tracks(trackIdx).xspeed = speed(2);
                obj.tracks(trackIdx).yspeed = speed(1);
                
                % replace predicted bounding box with detected
                % bounding box
                obj.tracks(trackIdx).bbox = bbox;
                obj.tracks(trackIdx).bboxProj = bboxProj;
                
                % update track's age
                obj.tracks(trackIdx).age = obj.tracks(trackIdx).age + 1;
                
                % update visibility
                obj.tracks(trackIdx).totalVisibleCount = ...
                    obj.tracks(trackIdx).totalVisibleCount + 1;
                
                % since we know that the track was visible for the previous
                % step, we set this variable to zero, as the name also suggests
                obj.tracks(trackIdx).consecutiveInvisibleCount = 0;
            end
        end
        
        function obj = updateUnassignedTracks(obj)
            % here we update the statistical parameters corresponding to
            % unassigned tracks
            for i = 1:length(obj.unassignedTracks) % for each unassigned track we loop
                ind = obj.unassignedTracks(i); % get the id of the track
                obj.tracks(ind).centroid = obj.tracks(ind).predCentroid;
                obj.tracks(ind).centroidObs = [-1,-1];
                obj.tracks(ind).centroidProj = obj.tracks(ind).predCentroidProj;
                obj.tracks(ind).centroidProjObs = [-1,-1];
                obj.tracks(ind).age = obj.tracks(ind).age + 1; % increment the age by one
                obj.tracks(ind).consecutiveInvisibleCount = ...
                    obj.tracks(ind).consecutiveInvisibleCount + 1; % since it is not observed, increase the invisible count
            end
        end
        
        function obj = deleteLostTracks(obj)
            if isempty(obj.tracks)
                return;
            end
            
            invisibleForTooLong = 5;
            ageThreshold = 3;
            
            % compute the fraction of the track's age for which it was visible
            ages = [obj.tracks(:).age];
            totalVisibleCounts = [obj.tracks(:).totalVisibleCount];
            visibility = totalVisibleCounts ./ ages;
            
            % find the indices of 'lost' tracks
            lostInds = (ages < ageThreshold & visibility < 0.1) | ...
                [obj.tracks(:).consecutiveInvisibleCount] >= invisibleForTooLong;
            
            % delete lost tracks
            %tracksLost = [tracksLost,tracks(lostInds)];
            obj.tracks = obj.tracks(~lostInds);
        end
        
        function obj = createNewTracks(obj)
            % here we create new tracks for the obbservations that are not
            % assigned to any existing tracks
            
            % first we filter out the necessary information corresponding to
            % unassigned observations. kind of information should be obvious
            % from looking at the name
            centroidsU = obj.centroids(obj.unassignedDetections, :);
            centroidsProjU = obj.centroidsProj(obj.unassignedDetections, :);
            bboxesU = obj.bboxes(obj.unassignedDetections, :);
            bboxesProjU = obj.bboxesProj(obj.unassignedDetections, :);
            areasU = obj.areas(obj.unassignedDetections, :);
            areasProjU = obj.areasProj(obj.unassignedDetections, :);
            contoursU = obj.contours(obj.unassignedDetections);
            contoursProjU = obj.contoursProj(obj.unassignedDetections);
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
                if(~isempty(obj.unassignedDetectionsToTracks) && sum(obj.unassignedDetectionsToTracks(:,2)==obj.unassignedDetections(i))>0)
                    trackIdx = obj.unassignedDetectionsToTracks(obj.unassignedDetectionsToTracks(:,2)==obj.unassignedDetections(i),1);
                    assignIdx = obj.assignments(obj.assignments(:,1)==trackIdx,2);
                    obj.tracks(trackIdx).kalmanFilter.State(obj.tracks(trackIdx).kalmanFilter.state_to_x) = obj.centroids(assignIdx,:);
                    
                    kalmanFilter = obj.tracks(trackIdx).kalmanFilter;
                    kalmanFilter.State(kalmanFilter.state_to_x) = centroid';
                    kalmanFilter = kalmanFilter.Update(centroid',double(area));
                    centroid = kalmanFilter.x';
                else
                    %kalmanFilter = SimpleKF([centroid,avgSpeed],0.6,[0.01,0.01],1e1,trackOption,double(area),1);
                    kalmanFilter = SimpleKF(centroid,0.6,[0.01,0.01],1e1,obj.trackOption,double(area),1);
                end
                
                % to track on the projected plane, we initialize another kalman
                % filter. The order of the input are the same as the previous
                % one.
                %kalmanFilterProj = SimpleKF(centroidProj,(obj.cbh^2/300^2)*1e1,(obj.cbh^2/300^2)*[1e3,1e2],(obj.cbh^2/300^2)*1e3,trackOption,double(areaProj),1);
                %kalmanFilterProj = SimpleKF(centroidProj,(obj.cbh^2/300^2)*1e2,(obj.cbh^2/300^2)*[9e4,1e3],(obj.cbh^2/300^2)*2e4,obj.trackOption,double(areaProj),1);
                %kalmanFilterProj = SimpleKF(centroidProj,1e2,[9e4,1e3],2e4,trackOption,double(areaProj),1);
                %kalmanFilterProj = SimpleKF(centroidProj,0.6,[0.01,0.01],1e1,trackOption,double(areaProj),1);
                if(~isempty(obj.unassignedDetectionsToTracks) && sum(obj.unassignedDetectionsToTracks(:,2)==obj.unassignedDetections(i))>0)
                    trackIdx = obj.unassignedDetectionsToTracks(obj.unassignedDetectionsToTracks(:,2)==obj.unassignedDetections(i),1);
                    assignIdx = obj.assignments(obj.assignments(:,1)==trackIdx,2);
                    obj.tracks(trackIdx).kalmanFilterProj.State(obj.tracks(trackIdx).kalmanFilterProj.state_to_x) = obj.centroidsProj(assignIdx,:);
                    
                    kalmanFilterProj = obj.tracks(trackIdx).kalmanFilterProj;
                    kalmanFilterProj.State(kalmanFilterProj.state_to_x) =  centroidProj';
                    kalmanFilterProj = kalmanFilterProj.Update(centroidProj',double(areaProj));
                    centroidProj = kalmanFilterProj.x';
                else
                    %kalmanFilterProj = SimpleKF([centroidProj,avgSpeedProj],(obj.cbh^2/300^2)*1e2,(obj.cbh^2/300^2)*[9e4,1e3],(obj.cbh^2/300^2)*2e4,trackOption,double(areaProj),1);
                    kalmanFilterProj = SimpleKF(centroidProj,(obj.cbh^2/300^2)*1e2,(obj.cbh^2/300^2)*[9e4,1e3],(obj.cbh^2/300^2)*2e4,obj.trackOption,double(areaProj),1);
                    %kalmanFilterProj = SimpleKF_v2(centroidProj,0.9,[0.01,0.15],1e1,trackOption,obj.cbh,1);
                end
                
                % create a new track
                newTrack = struct(...
                    'id', obj.nextId, ... % id of the track
                    'cbh', obj.cbh, ... % effective cbh assumed for the track
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
                    'obsSpeedProj', [0,0], ... % a variable keeping the 2D observed speed in projected plane
                    'area', double(area), ... % a variable keeping the area of the tracked object
                    'areaProj', double(areaProj), ... % a variable keeping the area of the tracked object
                    'contour', contour, ... % a variable keeping the relevant contour points on the image plane
                    'contourProj', contourProj, ... % a variable keeping the relevant contour points on the real world coordinates
                    'speedFilter', cloudTrackingCV.filterInit(10,2), ... % a filter for the observed speed on the image plane
                    'positionFilter', cloudTrackingCV.filterInit(10,2), ... % a filter for the position on the image plane
                    'centroidFilter', cloudTrackingCV.filterInit(1,2), ... % a filter for the centroid on the image plane
                    'centroidProjFilter', cloudTrackingCV.filterInit(1,2), ... % a filter for the centroid on the real plane
                    'massFilter', cloudTrackingCV.filterInit(10,1), ... % a filter on the predicted mass (not used)
                    'massFilterProj', cloudTrackingCV.filterInit(10,1), ... % a filter on the predicted mass (not used)
                    'centroid', centroid, ... % the representative position of the center of the object on the image plane
                    'centroidObs', centroid, ... % the representative observed position of the center of the object on the image plane
                    'centroidProj', centroidProj, ... % the representative position of the center of the object on the projected plane
                    'centroidProjObs', centroidProj, ... % the representative observed position of the center of the object on the projected plane
                    'predCentroid', centroid, ... % prediction on the representative position of the center of the object on the image plane
                    'predCentroidProj', centroidProj, ... % prediction on the representative position of the center of the object on the projected plane
                    'prevCentroid', centroid, ... % previous centroid position on the image plane
                    'prevCentroidProj', centroidProj); % previous centroid position on the projected plane
                
                
                % add it to the array of tracks
                obj.tracks(end + 1) = newTrack;
                
                % increment the next id
                obj.nextId = obj.nextId + 1;
            end
        end
        
        function obj = createNewPerfs(obj)
            persOcc = repmat(obj.occState,size(obj.occlusionValues));
            persIrr =  obj.irradianceValues_cs' ./ (persOcc*obj.attOcc + (1-persOcc)*obj.attCs);
            newPerf = struct(...
                'frameDateNum', datenum(obj.frameDateVec(1:6)'),...
                'frameCount', obj.frameCount,...
                'times', obj.irradianceTimes',...
                'persOcc', persOcc,...
                'algOcc', obj.occlusionValues,...
                'mesOcc', obj.occValues',...
                'persIrr', persIrr,...
                'algIrr', obj.irradianceValues_est',...
                'csIrr', obj.irradianceValues_cs',...
                'mesIrr', obj.irradianceValues');
            obj.perfs(end+1) = newPerf;
            %'persIrr', repmat(obj.irrMs,size(obj.occlusionValues)),...
            %'persOcc', repmat(obj.occState,size(obj.occlusionValues)),...
        end
        
        function displayTrackingResults(obj)
            %         % convert the frame and the mask to uint8 RGB
            %         frame = im2uint8(frame);
            %
            %         mask = uint8(repmat(mask, [1, 1, 3])) .* 255;
            %
            if(~obj.predict)
                minVisibleCount = 1;
                if ~isempty(obj.tracks)
                    
                    % noisy detections tend to result in short-lived tracks
                    % only display tracks that have been visible for more than
                    % a minimum number of frames.
                    reliableTrackInds = ...
                        [obj.tracks(:).totalVisibleCount] > minVisibleCount;
                    reliableTracks = obj.tracks(reliableTrackInds);
                    
                    % display the objects. If an object has not been detected
                    % in this frame, display its predicted bounding box.
                    if ~isempty(reliableTracks)
                        % get bounding boxes
                        obj.bboxes = cat(1, reliableTracks.bbox);
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
                            obj.mask = objectAnnotation(obj.mask, 'Rectangle', obj.bboxes(predictedTrackInds,:),bboxColor,3);
                            obj.mask = objectAnnotation(obj.mask, 'Text', [ids(predictedTrackInds,:),obj.bboxes(predictedTrackInds,1:2)],bboxColor,3);
                        end
                        if (sum(trackedTrackInds)>0)
                            bboxColor = 'y';
                            % draw on the mask
                            obj.mask = objectAnnotation(obj.mask, 'Rectangle', obj.bboxes(trackedTrackInds,:),bboxColor,3);
                            obj.mask = objectAnnotation(obj.mask, 'Text', [ids(trackedTrackInds,:),obj.bboxes(trackedTrackInds,1:2)],bboxColor,3);
                        end
                        if (sum(newTrackInds)>0)
                            bboxColor = 'g';
                            % draw on the mask
                            obj.mask = objectAnnotation(obj.mask, 'Rectangle', obj.bboxes(newTrackInds,:),bboxColor,3);
                            obj.mask = objectAnnotation(obj.mask, 'Text', [ids(newTrackInds,:),obj.bboxes(newTrackInds,1:2)],bboxColor,3);
                        end
                        
                        % draw on the mask
                        %mask = objectAnnotation(mask, 'Rectangle', bboxes,bboxColor,4);
                        
                        % draw on the frame
                        obj.frame = cv.drawContours(obj.frame,obj.contours,'Color',[255 255 0],'Thickness',3);
                        obj.predFrame = cv.drawContours(obj.predFrame,obj.contours,'Color',[0 255 0],'Thickness',3);
                    end
                end
                allPredContours = [obj.predCloudContours,obj.predShadowContours];
                colors = [repmat('b',size(obj.predCloudContours)),repmat('g',size(obj.predShadowContours))];
                % draw on the mask
                
                %mask = objectAnnotation(mask, 'Rectangle', bboxes);
                
                % draw on the frame
                %frame = cv.drawContours(frame,contours,'Color',[255 255 0]);
                
                % display the mask and the frame
                obj.maskPlayer.Step(obj.mask);
                obj.videoPlayer.Step(obj.frame);
                obj.predPlayer.Step(obj.predFrame);
                obj.mapPlayer = obj.mapPlayer.Step(allPredContours,colors);
                if(obj.gifMode)
                    frame1 = obj.videoPlayer.GetFrame();
                    frame2 = obj.maskPlayer.GetFrame();
                    frame3 = obj.predPlayer.GetFrame();
                    %pause(0.02);
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
                minVisibleCount = 1;
                if ~isempty(obj.tracks)
                    
                    % noisy detections tend to result in short-lived tracks
                    % only display tracks that have been visible for more than
                    % a minimum number of frames.
                    reliableTrackInds = ...
                        [obj.tracks(:).totalVisibleCount] > minVisibleCount;
                    reliableTracks = obj.tracks(reliableTrackInds);
                    
                    % display the objects. If an object has not been detected
                    % in this frame, display its predicted bounding box.
                    if ~isempty(reliableTracks)
                        % get bounding boxes
                        obj.bboxes = cat(1, reliableTracks.bbox);
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
                            obj.mask = objectAnnotation(obj.mask, 'Rectangle', obj.bboxes(predictedTrackInds,:),bboxColor,3);
                            %obj.mask = objectAnnotation(obj.mask, 'Text', [ids(predictedTrackInds,:),obj.bboxes(predictedTrackInds,1:2)],bboxColor,3);
                        end
                        if (sum(trackedTrackInds)>0)
                            bboxColor = 'y';
                            % draw on the mask
                            obj.mask = objectAnnotation(obj.mask, 'Rectangle', obj.bboxes(trackedTrackInds,:),bboxColor,3);
                            %obj.mask = objectAnnotation(obj.mask, 'Text', [ids(trackedTrackInds,:),obj.bboxes(trackedTrackInds,1:2)],bboxColor,3);
                        end
                        if (sum(newTrackInds)>0)
                            bboxColor = 'g';
                            % draw on the mask
                            obj.mask = objectAnnotation(obj.mask, 'Rectangle', obj.bboxes(newTrackInds,:),bboxColor,3);
                            %obj.mask = objectAnnotation(obj.mask, 'Text', [ids(newTrackInds,:),obj.bboxes(newTrackInds,1:2)],bboxColor,3);
                        end
                        
                        % draw on the mask
                        %obj.mask = objectAnnotation(obj.mask, 'Rectangle', obj.bboxes,bboxColor,4);
                        
                        % draw on the frame
                        obj.frame = cv.drawContours(obj.frame,obj.contours,'Color',[255 255 0],'Thickness',3);
                        obj.predFrame = cv.drawContours(obj.predFrame,obj.contours,'Color',[0 255 0],'Thickness',3);
                    end
                end
                %allPredContours = [predCloudContours,predShadowContours];
                %colors = [repmat('b',size(predCloudContours)),repmat('g',size(predShadowContours))];
                % draw on the mask
                
                %mask = objectAnnotation(mask, 'Rectangle', bboxes);
                
                % draw on the frame
                %frame = cv.drawContours(frame,contours,'Color',[255 255 0]);
                
                % display the mask and the frame
                obj.maskPlayer.Step(obj.mask,obj.frameDateStr);
                %obj.videoPlayer.Step(frame);
                obj.predPlayer.Step(obj.predFrame,obj.frameDateStr);
                %obj.mapPlayer = obj.mapPlayer.Step(allPredContours,colors);
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
                % if we are dumping to a webseite.
                textColor = 'g';
                irradianceVal = obj.cbh;
                %irradianceValues = obj.irradianceValues(end);
                value_positions = [repmat(70,1,length(irradianceVal));70:35:70+((length(irradianceVal)-1)*35)];
                obj.frame = objectAnnotation(obj.frame, 'Text', [irradianceVal, value_positions'], textColor,3);
                obj.videoPlayer.Step(obj.frame,obj.frameDateStr);
                obj.plotPlayerGHI.Step([obj.irradianceTimes,obj.irradianceValues,obj.irradianceTimes_cs,obj.irradianceValues_cs],strcat('Time Relative to ', obj.frameDateStr, ' (Min.)') ,'Irradiance (W/m2)','Irradiance Value Changes',['Measured ';'Clear-Sky']);
                obj.plotPlayerPV.Step([obj.irradianceTimes,obj.irradianceValues,obj.irradianceTimes,obj.irradianceValues_est],strcat('Time Relative to ', obj.frameDateStr, ' (Min.)') ,'Irradiance (W/m2)','Irradiance Value Changes',['Measured ';'Predicted']);
                %obj.plotPlayerPV.Step([obj.powerTimes,obj.powerValues],strcat('Time Relative to ', obj.frameDateStr, ' (Min.)') ,'Power (W)','Power Production Value Changes');
                %obj.plotPlayerTemp.Step([obj.tempTimes,obj.tempValues],strcat('Time Relative to ', obj.frameDateStr, ' (Min.)') ,'Temperature (degrees)','Temperature Value Changes');
                obj.plotPlayerTemp.Step([obj.occTimes,obj.occValues,obj.occlusionTimes',obj.occlusionValues'],strcat('Time Relative to ', obj.frameDateStr, ' (Min.)') ,'Occlusion (binary)','Occlusion Value Changes',['Measured ';'Predicted']);
                if(obj.vid_r)
                    ff11 = obj.videoPlayer.GetFrame();
                    ff12 = obj.maskPlayer.GetFrame();
                    ff13 = obj.predPlayer.GetFrame();
                    ff21 = obj.plotPlayerGHI.GetFrame();
                    ff22 = obj.plotPlayerPV.GetFrame();
                    ff23 = obj.plotPlayerTemp.GetFrame();
                    ff = [ff11,ff12,ff13;ff21,ff22,ff23];
                    writeVideo(obj.vidObj, ff);
                end
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
        
        function obj = closeVideoObj(obj)
            if(obj.vid_r)
                close(obj.vidObj);
            end
        end
        
        function tracks = getTracks(obj)
            tracks = obj.tracks;
        end
        
        function perfs = getPerfs(obj)
            perfs = obj.perfs;
        end
        
        function reader = getReader(obj)
            reader = obj.reader;
        end
        
        function frameCount = getFrameCount(obj)
            frameCount = obj.frameCount;
        end
        
        function finalFrame = getFinalFrame(obj)
            finalFrame = obj.finalFrame;
        end
        
        function sunEliminationCycle = getSunEliminationCycle(obj)
            sunEliminationCycle = obj.sunEliminationCycle;
        end
        
        function frame = getFrame(obj)
            frame = obj.frame;
        end
        
        function frameRate = getFrameRate(obj)
            frameRate = obj.frameRate;
        end
        
        function sunPosition = getSunPlacement(obj)
            sunPosition = obj.sunPosition;
        end
    end
    
    methods(Static)
        function tracks = initializeTracks()
            % create an empty array of tracks
            tracks = struct(...
                'id', {}, ...   % id of the track
                'cbh', {}, ...  % effective cbh assumed for the track
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
                'obsSpeedProj', {}, ... % a variable keeping the 2D observed speed in projected plane
                'area', {}, ... % a variable keeping the area of the tracked object
                'areaProj', {}, ... % a variable keeping the area of the tracked object
                'contour', {}, ... % a variable keeping the relevant contour points on the image plane
                'contourProj', {}, ... % a variable keeping the relevant contour points on the real world coordinates
                'speedFilter', {}, ... % a filter for the observed speed on the image plane
                'positionFilter', {}, ... % a filter for the position on the image plane
                'centroidFilter', {}, ... % a filter for the centroid on the image plane
                'centroidProjFilter', {}, ... % a filter for the centroid on the real plane
                'massFilter', {}, ... % a filter on the predicted mass (not used)
                'massFilterProj', {}, ... % a filter on the predicted mass (not used)
                'centroid', {}, ... % the representative position of the center of the object on the image plane
                'centroidObs', {}, ... % the representative observed position of the center of the object on the image plane
                'centroidProj', {}, ... % the representative position of the center of the object on the projected plane
                'centroidProjObs', {}, ... % the representative observed position of the center of the object on the projected plane
                'predCentroid', {}, ... % prediction on the representative position of the center of the object on the image plane
                'predCentroidProj', {}, ... % prediction on the representative position of the center of the object on the projected plane
                'prevCentroid', {}, ... % previous centroid position on the image plane
                'prevCentroidProj', {}); % previous centroid position on the projected plane
            
        end
        
        function perfs = initializePerfs()
            perfs = struct(...
                'frameDateNum', {},...
                'frameCount', {},...
                'times', {},...
                'persOcc', {},...
                'algOcc', {},...
                'mesOcc', {},...
                'persIrr', {},...
                'algIrr', {},...
                'csIrr', {},...
                'mesIrr', {});
        end
        
        function pointReal = getImage2Real_v2(m,calibModel,ecalibModel,cbh)
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
        
        function pointImage = getReal2Image_v2(m,calibModel,ecalibModel,cbh)
            % here we project a set of points given in m back onto the image plane. The used convention in opencv is [column,row] kind of
            % represantation for pixel points, here however we change that to [row,column]
            if(size(ecalibModel.R,1)==3)
                pointTemp = world2cam(ecalibModel.Rinv*[m';-cbh*ones(1,size(m,1))],calibModel)';
            else
                pointTemp = world2cam([ecalibModel.R'*m';-cbh*ones(1,size(m,1))],calibModel)'; % we get the corresponding direction on a unit-sphere at first.
            end
            pointImage = [pointTemp(:,2),pointTemp(:,1)]-1;
        end
        
        function pointReal = getImage2Real(m,calibModel,ecalibModel,cbh)
            % here we project a set of points given in m onto a plane that is cbh meters away from the camera the used convention in opencv is [column,row] kind of
            % represantation for pixel points, here however we change that to [row,column]
            R_earth = 6378100;
            pointTemp = cam2world([m(:,2),m(:,1)]',calibModel); % we get the corresponding direction on a unit-sphere at first.
            pointReal = zeros(size(pointTemp,2),2); % we crate the array that will contain the projected positions
            if(size(ecalibModel.R,1)==3)
                pointTemp = ecalibModel.R*pointTemp;
                alphas = (pi/2) - atan2(-pointTemp(3,:),sqrt(pointTemp(1,:).^2 + pointTemp(2,:).^2));
                heights = (-2*R_earth*cos(alphas) + sqrt((2*R_earth*cos(alphas)).^2 + 4*(cbh^2+2*cbh*R_earth)))/2;
                pointReal = (pointTemp(1:2,:).*repmat(heights,2,1))';
            else
                pointReal = ((ecalibModel.R*pointTemp(1:2,:)).*repmat((-cbh./pointTemp(3,:)),2,1))'; % here we map the point on the unit-sphere onto the plane
            end
        end
        
        function pointImage = getReal2Image(m,calibModel,ecalibModel,cbh)
            % here we project a set of points given in m back onto the image plane. The used convention in opencv is [column,row] kind of
            % represantation for pixel points, here however we change that to [row,column]
            R_earth = 6378100;
            if(size(ecalibModel.R,1)==3)
                heights = sqrt((R_earth+cbh)^2-(m(:,1).^2 + m(:,2).^2))-R_earth;
                pointTemp = world2cam(ecalibModel.Rinv*[m';-heights'],calibModel)';
            else
                pointTemp = world2cam([ecalibModel.R'*m';-cbh*ones(1,size(m,1))],calibModel)'; % we get the corresponding direction on a unit-sphere at first.
            end
            pointImage = [pointTemp(:,2),pointTemp(:,1)]-1;
        end
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
    end
end