%% Cloud tracking using kalman filter.
%
% Multiple object tracking with kalman filter for could tracking:
%Inputs:
% video_name: path to the input video file
% record: option indicating that a record should be taken or not
% verbose: option indication that the video should be seen or not
% finalFrame: number of frames of the video, that the algorithm will run
% on.
%Outputs:
% recording_obs: Structure array that keeps the set of all observations. 
% Kept mostly for evaluation and plotting purposes.
% recording_pred: Structure array that keeps the set of all predictions for
% all the timesteps. Kept mostly for evaluation and plotting purposes.
% db_objs: Recordingos of the observations.
% db_ass: Recordings of the observations to previously tracked objects.
% B. Zeydan, 27. Feb. 2013

function [recording_obs,recording_pred,db_objs,db_ass] = multiObjectTrackingCV(video_name,record,verbose,finalFrame)
% create system objects used for 2reading video, detecting moving objects,
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
obj = setupSystemObjects(video_name,verbose);

predictionCycle = 300;
updateCycle = 1;

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
% detect moving objects, and track them across video frames
while obj.frameCount < finalFrame %~isDone(obj.reader)
    frame = readFrame();
    camParams = camParamsInit();
    
    if((mod(obj.frameCount,predictionCycle) < updateCycle) || (obj.frameCount < 3))
        [contours,areas,centroids, bboxes, mask] = detectObjects(frame);
        
        db_objs{obj.frameCount} = struct('contours',contours,...
            'areas',areas,...
            'centroids',centroids,...
            'bboxes',bboxes);
        
        predictNewLocationsOfTracks();
        if(record && (predictionCycle >= updateCycle))
            recording_pred{obj.frameCount} = tracks;
        end
        [assignments, unassignedTracks, unassignedDetections] = ...
            detectionToTrackAssignment();
        
        db_ass{obj.frameCount} = struct('assignments', assignments,...
            'unassignedTracks', unassignedTracks,...
            'unassignedDetections', unassignedDetections);
        
        updateAssignedTracks();
        updateUnassignedTracks();
        deleteLostTracks();
        createNewTracks();
    else
        mask = getBinaryImage(frame,maskMethod);
        mask = medfilt2(mask,[15 15]);
        
        %mask = imclose(mask, strel('square', 30));
        %mask = imfill(mask, 'holes');
        predictNewLocationsOfTracks();
        updateBlindTracks();
        deleteLostTracks();
    end
    
    % for recording purposes...
    if(record)
        if(obj.frameCount == 1)
            recording_pred{obj.frameCount} = tracks;
        end
        recording_obs{obj.frameCount} = tracks;
    end
    if(verbose)
        displayTrackingResults();
    end
end

    function obj = setupSystemObjects(video_name,verbose)
        % Initialize Video I/O
        % Create objects for reading a video from a file, drawing the tracked
        % objects in each frame, and playing the video.
        
        % create a video file reader
        obj.reader = VideoReader(video_name);
        obj.frameCount = 0;
        
        % create two video players, one to display the video,
        % and one to display the foreground mask
        if(verbose)
            obj.videoPlayer = VidPlayer([20, 400, 700, 400]);
            obj.maskPlayer = VidPlayer([740, 400, 700, 400]);
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

    function mask = getBinaryImage(frame, type)
        
        r_b = double(frame(:,:,1))./double(frame(:,:,3)); % R/B thresholding applied
        mask = r_b > 0.72;
        %mask = im2bw(frame,0.45);
        
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
                mask = objectAnnotation(mask, 'Rectangle', bboxes,bboxColor);
                
                % draw on the frame
                frame = cv.drawContours(frame,contours,'Color',[255 255 0]);
                
            end
        end
        
        % draw on the mask
        
        %mask = objectAnnotation(mask, 'Rectangle', bboxes);
        
        % draw on the frame
        %frame = cv.drawContours(frame,contours,'Color',[255 255 0]);
        
        % display the mask and the frame
        obj.maskPlayer.Step(mask);
        obj.videoPlayer.Step(frame);
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
        
        cbh = 5000;
        
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
