function [recording_tracks,perfs] = multiObjectTrackingCV_class(vid_n,vid_r,image_folder,record,verbose,calibModel,ecalibModel,model3D,imageMask,sunPattern,modelImage,live,track,timeSeriesGHI,timeSeriesGHI_cs,timeSeriesPV,timeSeriesTemp,timeSeriesOcc,timeSeriesCBH)
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
recording_tracks = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the variables that determine the options and parameters applied on the
% tracking algorithm.

cbh = 1500; % the assumed cloud base height parameter that determines the projection performance of the algorithm.
ModPlayerOption = 'map'; % the representation of clouds on the rea-sized model is done using this opttion. rectangle is another possibility.
maskMethod = 'body'; % identification type of clouds is done by this variable. either they are identified by their body or the contours (edges)
trackOption = 'speed'; % the option that is applied for the kalman filter, for other possibilities please look inside the object implementation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% here we set-up the relevant objects.
cloudTracking = cloudTrackingCV(vid_n,vid_r,image_folder,live,verbose,calibModel,ecalibModel,model3D,imageMask,modelImage,cbh,ModPlayerOption,track,timeSeriesGHI,timeSeriesGHI_cs,timeSeriesPV,timeSeriesTemp,timeSeriesOcc,timeSeriesCBH,maskMethod,trackOption,sunPattern);
% in case we go with the validation mode and we are recording the
% performance of the algorithm, the pstate variables of all the steps are
% recorded in these variables
%obj.finalFrame = 200;
if(isempty(cloudTracking.getReader()))
    disp('Camera does not respond. Terminating...');
    return;
end

avg_t = 0;

% This is the main loop of the program
while cloudTracking.getFrameCount() < cloudTracking.getFinalFrame()
    tst = tic;
    cloudTracking = cloudTracking.readFrame();  % readFrame function returns the image, that is in order, from the buffer folder if this is not live, or the recording folder that is determined during the setup
    cloudTracking = cloudTracking.prepareData();
    
    % This takes time, because it matches the sun pattern in the image, so
    % we run it not so frequently...
    cloudTracking = cloudTracking.getSunPositionReal();
    cloudTracking = cloudTracking.getSunPatch();
    %if((sum(cloudTracking.getSunPlacement())<0) || mod(cloudTracking.getFrameCount(),cloudTracking.getSunEliminationCycle())>=(cloudTracking.getSunEliminationCycle()-cloudTracking.getFrameRate()))
    %    cloudTracking = cloudTracking.getSunPosition();
    %    cloudTracking = cloudTracking.getSunPatch();
    %    cloudTracking = cloudTracking.getShadowDisplacement();
    %end
    
    if(isempty(cloudTracking.getFrame()))    % if the folder does not exist or the frame is not correct, we immediately terminate
        disp('Camera does not respond. Terminating...');
        close all;
        return;
    end
    
    cloudTracking = cloudTracking.eliminateSun();
    cloudTracking = cloudTracking.getBinaryImage('body');
    cloudTracking = cloudTracking.detectObjects(); % this function takes in the RGB frame and :
    % - runs the cloud decision algorithm to get a binary image
    % - runs a clustering algorithm followed by a contour extraction
    % - runs simple heuristic methods on the contour pixels to extract
    %   - the representative points for the model based tracking
    %   - the projected representative point for physical interpretation
    %   - bounding boxes for the objects found
    % - returns all the results along with the binary image to be used for
    % displaying
    
    
    cloudTracking = cloudTracking.predictNewLocationsOfTracks();    % in this function we predict the one-step future of the existing model calling the predict method of
    % existing Kalman filters.
    
    cloudTracking = cloudTracking.predictNewLocationsOfTracksIntoFuture();
    
    if(verbose && ~record)
        %[predTracks,obj.occlusionValues,predCloudContours,predShadowContours] = predictNewLocationsOfTracksIntoFuture(obj.timeHorizon,obj.occlusionValues);
        cloudTracking = cloudTracking.getPredictedImage();
        %if(sum(obj.occlusionValues)>0)
        %    obj.occlusionValues
        %end
        %[predCloudContours,predShadowContours] = getPredictedScenario(predTracks);
    end
    
    cloudTracking = ...
        cloudTracking.detectionToTrackAssignment(); % for every frame that is processed, there is a set of detections. Also for the frames previous to the current one, there are the tracks.
    % This function establishes the association between the observations and tracks.
    
    cloudTracking = cloudTracking.updateAssignedTracks(); % this is the function that updates the tracks according to assignments, it also updates some statistical variables related to update rate
    cloudTracking = cloudTracking.updateUnassignedTracks(); % this is the funstion that updates tracks corresponding to unassigned tracks, it just updates statistical variables related to update rate
    cloudTracking = cloudTracking.deleteLostTracks(); % here we delete tracks that remained unassigned for some number of steps
    cloudTracking = cloudTracking.createNewTracks(); % here we create new tracks corresponding to unassigned observations
    cloudTracking = cloudTracking.createNewPerfs();
    %     if(~record)
    %         occluded = getOcclusion(predTracks,model3D);
    %         if(occluded)
    %             occluded
    %         end
    %     end
    cloudTracking = cloudTracking.createNewPerfs();
    
    if(record)
        if(isempty(recording_tracks))
            recording_tracks{1} = cloudTracking.getTracks();
        else
            recording_tracks{end+1} = cloudTracking.getTracks();  % this is the database where we record the predictions of the model before the updates.
        end
    end
    % for recording purposes...
    if(verbose && ~record)
        cloudTracking.displayTrackingResults(); % here we display the status and results
    end
    
    if(avg_t == 0)
        avg_t = toc(tst);
    else
        avg_t = 0.8*avg_t + 0.2*toc(tst);
    end
    
    if(mod(cloudTracking.getFrameCount(),100)==5 && ~verbose)
        %obj.frameCount
        if(avg_t~=0)
            disp(strcat('Remaining Time: ~',num2str(round(((cloudTracking.getFinalFrame()-cloudTracking.getFrameCount())/cloudTracking.getFrameRate())*avg_t/60)),' minutes.'))
        end
    elseif(mod(cloudTracking.getFrameCount(),100)==5 && verbose)
        if(avg_t~=0)
            disp(strcat('Average time per frame: ~',num2str(avg_t),{' seconds.'}));
        end
    end
end
perfs = cloudTracking.getPerfs();
cloudTracking.closeVideoObj();
end