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


function [recording_obs,recording_pred] = multiObjectTrackingCV_Correlation(video_name,record,verbose,finalFrame)
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
obj = setupSystemObjects(video_name);

inf = info(obj.reader);
frameSize = inf.VideoSize;
R = 100;
y_edges = [1:R:frameSize(1)-R+1];
x_edges = [1:R:frameSize(2)-R+1];
[bbx,bby] = meshgrid(y_edges,x_edges);
size([bbx,bby])
mesh_lut = [bbx(:),bby(:),(R-1)*ones(size([bbx(:),bby(:)]))];
mesh_cent_lut = [bbx(:),bby(:)] + 0.5*R*ones(size([bbx(:),bby(:)]));
tracks = initializeTracks(); % create an empty array of tracks
mild_noise = 0.0001*rand(R,R);
nextId = 1; % ID of the next track

trackOption = 'pose_v';
areas = [];
centroids = [];
bboxes = [];
mask = [];
obs = [];
if(record)
    tracks2 = initializeTracks(); % create an empty array of tracks
    
    nextId2 = 1; % ID of the next track
    
    trackOption2 = 'pose_v';
    areas2 = [];
    centroids2 = [];
    bboxes2 = [];
    mask2 = [];
end
frameCount = 1;
% detect moving objects, and track them across video frames
while frameCount < finalFrame %~isDone(obj.reader)
    frame = readFrame();
    camParams = camParamsInit();
    
    
    mask = getBinaryImage(frame,'body');
    %mask = medfilt2(mask,[5 5]);
    
    if(frameCount == 1)
        prev_mask = mask;
    end
    
    %tolmask = zeros(size(mask)+6*[R,R]);
    %tolmask(3*R:3*R+size(mask,1),3*R:3*R+size(mask,2)) = mask;
    %[max_val, max_ind] = cellfun(@(x) max(abs(normxcorr2(prev_mask(x(2):x(2)+x(4),x(1):x(1)+x(3))+mild_noise,tolmask(3*R+x(2)-3*R:3*R+x(2)+x(4)+3*R,3*R+x(1)-3*R:3*R+x(1)+x(3)+3*R)))) ,mat2cell(mesh_lut,ones(size(mesh_lut,1),1)), 'UniformOutput',false);
    [max_val, max_ind] = cellfun(@(x) max(abs(normxcorr2(prev_mask(x(2):x(2)+x(4),x(1):x(1)+x(3))+mild_noise,mask))) ,mat2cell(mesh_lut,ones(size(mesh_lut,1),1)), 'UniformOutput',false);
    [max_ind_x,max_ind_y] = cellfun(@(x) ind2sub(frameSize,x(1)),max_ind,'UniformOutput',false);
    
    ll = [max_ind_x{:};max_ind_y{:}];
    
    arrows = [mesh_cent_lut,ll'];
    
    if(verbose)
        displayTrackingResults();
    end
    frameCount = frameCount + 1;
end

    function frame = negativeFrame(frame)
        i_frame = ones(size(frame));
        frame = i_frame-frame;
    end

    function obj = setupSystemObjects(video_name)
        % Initialize Video I/O
        % Create objects for reading a video from a file, drawing the tracked
        % objects in each frame, and playing the video.
        
        % create a video file reader
        obj.reader = vision.VideoFileReader(video_name);
        
        % create two video players, one to display the video,
        % and one to display the foreground mask
        obj.videoPlayer = vision.VideoPlayer('Position', [20, 400, 700, 400]);
        obj.maskPlayer = vision.VideoPlayer('Position', [740, 400, 700, 400]);
        obj.videoShapeInserter = vision.ShapeInserter('Shape','Lines','BorderColor','Custom','CustomBorderColor',[255,255,0]);
        obj.maskShapeInserter = vision.ShapeInserter('Shape','Lines','BorderColor','Custom','CustomBorderColor',[255,255,0]);
        obj.videoMarkerInserter = vision.MarkerInserter('Shape','Star','BorderColor','Custom','CustomBorderColor',[255,0,0]);
        obj.maskMarkerInserter = vision.MarkerInserter('Shape','Star','BorderColor','Custom','CustomBorderColor',[255,0,0]);
        
        % Create system objects for foreground detection and blob analysis
        
        % The foreground detector is used to segment moving objects from
        % the background. It outputs a binary mask, where the pixel value
        % of 1 corresponds to the foreground and the value of 0 corresponds
        % to the background.
        
        obj.detector = vision.ForegroundDetector('NumGaussians', 3, ...
            'NumTrainingFrames', 5, 'MinimumBackgroundRatio', 0.001);
        
        % Connected groups of foreground pixels are likely to correspond to moving
        % objects.  The blob analysis system object is used to find such groups
        % (called 'blobs' or 'connected components'), and compute their
        % characteristics, such as area, centroid, and the bounding box.
        
        obj.blobAnalyser = vision.BlobAnalysis('BoundingBoxOutputPort', true, ...
            'AreaOutputPort', true, 'CentroidOutputPort', true, ...
            'MinimumBlobArea', 900);
    end

    function tracks = initializeTracks()
        % create an empty array of tracks
        tracks = struct(...
            'id', {}, ...
            'bbox', {}, ...
            'kalmanFilter', {}, ...
            'simpleKF', {}, ...
            'age', {}, ...
            'totalVisibleCount', {}, ...
            'consecutiveInvisibleCount', {},...
            'xspeed', {}, ...
            'yspeed', {}, ...
            'xrspeed', {}, ...
            'yrspeed', {}, ...
            'area', {}, ...
            'speedFilter', {}, ...
            'positionFilter', {}, ...
            'massFilter', {}, ...
            'centroid', {});
        
    end
    
    function frame = readFrame()
        frame = obj.reader.step();
    end

    function mask = getBinaryImage(frame, type)
        
        %r_b = frame(:,:,1)./frame(:,:,3); % R/B thresholding applied
        %mask = r_b > 0.72;
        mask = frame(:,:,1);
        %mask = im2bw(frame,0.45);
        
        cloudLevel = sum(sum(mask))/numel(mask);
        
        % if high cloud density, than we follow the sky segments..
        if(cloudLevel > 0.8)
            mask = r_b < 0.72;
            %mask = ones(size(mask)) - mask;
        end
        
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
            predictedCentroid = predict(tracks(i).kalmanFilter);
            
            tracks(i).simpleKF = tracks(i).simpleKF.Estimate();
            % shift the bounding box so that its center is at
            % the predicted location
            predictedCentroid = int32(predictedCentroid(1:2)) - bbox(3:4) / 2;
            
            tracks(i).bbox = [predictedCentroid(1:2), bbox(3:4)];
        end
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
            detectionToTrackAssignment()
        nTracks = length(tracks);
        nDetections = size(centroids, 1);
        
        % compute the cost of assigning each detection to each track
        cost = zeros(nTracks, nDetections);
        for i = 1:nTracks
            if(strcmp(trackOption,'areas'))
                cost(i, :) = distance(tracks(i).kalmanFilter, [centroids,areas]);
            else
                cost(i, :) = distance(tracks(i).kalmanFilter, centroids);
            end
        end
        
        % solve the assignment problem
        costOfNonAssignment = 10000;
        [assignments, unassignedTracks, unassignedDetections] = ...
            assignDetectionsToTracks(cost, costOfNonAssignment);
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
            tracks(trackIdx).centroid = [tracks(trackIdx).kalmanFilter.State(1),tracks(trackIdx).kalmanFilter.State(3)];
        end
        for i = 1:length(unassignedTracks)
            ind = unassignedTracks(i);
            if(ind > length(tracks))
                break;
            end
            tracks(ind).age = tracks(ind).age + 1;
            tracks(ind).consecutiveInvisibleCount = ...
                tracks(ind).consecutiveInvisibleCount + 1;
            tracks(ind).centroid = [tracks(ind).kalmanFilter.State(1),tracks(ind).kalmanFilter.State(3)];
        end
    end

    function updateAssignedTracks()
        numAssignedTracks = size(assignments, 1);
        frameSize = size(frame);
        for i = 1:numAssignedTracks
            trackIdx = assignments(i, 1);
            detectionIdx = assignments(i, 2);
            centroid = centroids(detectionIdx, :);
            bbox = bboxes(detectionIdx, :);
            area = areas(detectionIdx, :);
            
            %prevCentroid = tracks(trackIdx).kalmanFilter.State(1:2:3);
            % correct the estimate of the object's location
            % using the new detection
            
            pos_t = (centroid - (frameSize(1:2)/2));%*camParams.conversionFactor;
            
            %pos_t = centroid;
            tracks(trackIdx).positionFilter = filterAdd(pos_t,tracks(trackIdx).positionFilter);
            pos_t = tracks(trackIdx).positionFilter.value;
            
            tracks(trackIdx).simpleKF = tracks(trackIdx).simpleKF.Update(pos_t(1),pos_t(2),double(area));
            
            if(strcmp(trackOption,'areas'))
                correct(tracks(trackIdx).kalmanFilter, [centroid,area]);
            else
                correct(tracks(trackIdx).kalmanFilter, centroid);
            end
            
            % getting the speed estimate from the kalman...
            %speed = [tracks(trackIdx).simpleKF.x(2), tracks(trackIdx).simpleKF.x(5)];
            if(strcmp(trackOption,'pose_v'))
                speed = [tracks(trackIdx).kalmanFilter.State(2), tracks(trackIdx).kalmanFilter.State(4)];
            elseif(strcmp(trackOption,'pose_acc'))
                speed = [tracks(trackIdx).kalmanFilter.State(2), tracks(trackIdx).kalmanFilter.State(5)];
            end
            
            tracks(trackIdx).speedFilter = filterAdd(speed, tracks(trackIdx).speedFilter);
            speed = tracks(trackIdx).speedFilter.value;% * camParams.conversionFactor;
            
            tracks(trackIdx).xspeed = speed(2);
            tracks(trackIdx).yspeed = speed(1);
            
            tracks(trackIdx).centroid = centroid;
            
            tracks(trackIdx).massFilter = filterAdd(area,tracks(trackIdx).massFilter);
            tracks(trackIdx).area = tracks(trackIdx).massFilter.value;
            
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

    function updateUnassignedTracks()
        for i = 1:length(unassignedTracks)
            ind = unassignedTracks(i);
            tracks(ind).age = tracks(ind).age + 1;
            tracks(ind).consecutiveInvisibleCount = ...
                tracks(ind).consecutiveInvisibleCount + 1;
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

    function deleteLostTracks()
        if isempty(tracks)
            return;
        end
        
        invisibleForTooLong = 900;
        ageThreshold = 8;
        
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

    function createNewTracks()
        centroids = centroids(unassignedDetections, :);
        bboxes = bboxes(unassignedDetections, :);
        areas = areas(unassignedDetections, :);
        
        for i = 1:size(centroids, 1)
            
            centroid = centroids(i,:);
            bbox = bboxes(i, :);
            
            area = areas(i);
            
            % create a Kalman filter object
            %kalmanFilter = configureKalmanFilter('ConstantVelocity', ...
            %    centroid, [10, 25], [100, 25], 100);
            if(strcmp(trackOption,'areas'))
                state = [centroid(1), 0 , centroid(2), 0, area, 0];
                state = [centroid(1), 0 , centroid(2), 0, 0, 0];
                kf_dims = 6;
                P = 1*eye(kf_dims);
                Q = 0.1 * eye(kf_dims);
                Q(1,1) = 3;
                Q(3,3) = 3;
                Q(5,5) = 1;
                R = 20 * [1 0 0; 0 1 0; 0 0 1]; % we have 3 measurements, so this should be in this shape.
                A = [1 1 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 1 1 0 0; 0 0 0 0 1 1; 0 0 0 0 0 1]; %
                H = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
                kalmanFilter = vision.KalmanFilter('StateTransitionModel', A, 'MeasurementModel', H, ...
                    'State', state, 'StateCovariance', P, 'ProcessNoise', Q, 'MeasurementNoise', R);
            elseif(strcmp(trackOption,'pose_v'))
                kalmanFilter = configureKalmanFilter('ConstantVelocity', ...
                    centroid, [3, 3], [3, 3], 100);
            elseif(strcmp(trackOption,'pose_acc'))
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
                'xrspeed', 0, ...
                'yrspeed', 0, ...
                'area', double(area), ...
                'speedFilter', filterInit(10,2), ...
                'positionFilter', filterInit(10,2), ...
                'massFilter', filterInit(10,1), ...
                'centroid', centroid);
            
            %newTrack.rpixel = newTrack.xsensor/frameSize(2);
            %newTrack.rsensor = sqrt(newTrack.xsensor^2 + newTrack.ysensor^2);
            %newTrack.conversionFactor = sqrt((0.002/tan(45))^2 + (0.002^2 + 0.0013^2))/5000;
            % add it to the array of tracks
            tracks(end + 1) = newTrack;
            
            % increment the next id
            nextId = nextId + 1;
        end
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
                'id', nextId2, ...
                'bbox', bbox, ...
                'kalmanFilter', kalmanFilter, ...
                'simpleKF', simpleKF, ...
                'age', 1, ...
                'totalVisibleCount', 1, ...
                'consecutiveInvisibleCount', 0, ...
                'xspeed', 0, ...
                'yspeed', 0, ...
                'xrspeed', 0, ...
                'yrspeed', 0, ...
                'area', double(area), ...
                'speedFilter', filterInit(10,2), ...
                'positionFilter', filterInit(10,2), ...
                'massFilter', filterInit(10,1), ...
                'centroid', centroid);
            
            %newTrack.rpixel = newTrack.xsensor/frameSize(2);
            %newTrack.rsensor = sqrt(newTrack.xsensor^2 + newTrack.ysensor^2);
            %newTrack.conversionFactor = sqrt((0.002/tan(45))^2 + (0.002^2 + 0.0013^2))/5000;
            % add it to the array of tracks
            tracks2(end + 1) = newTrack;
            
            % increment the next id
            nextId2 = nextId2 + 1;
        end
    end

    function displayTrackingResults()
        % convert the frame and the mask to uint8 RGB
        frame = im2uint8(frame);
        
        mask = uint8(repmat(mask, [1, 1, 3])) .* 255;
        
        
        frame = step(obj.videoShapeInserter, frame, uint32(arrows));
        frame = step(obj.videoMarkerInserter, frame, uint32(mesh_cent_lut));
        % draw on the mask
        mask = step(obj.maskShapeInserter, mask, uint32(arrows));
        mask = step(obj.maskMarkerInserter, mask, uint32(mesh_cent_lut));
        
        % display the mask and the frame
        obj.maskPlayer.step(mask);
        obj.videoPlayer.step(frame);
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
