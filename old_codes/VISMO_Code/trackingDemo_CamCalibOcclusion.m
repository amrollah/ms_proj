%% Demo of cloud tracking code using camera inputs and calibration.
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
% B. Zeydan, 1. May. 2013

%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/05/14'; % the name of the stream file (video or a collection of images)
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/06/28';
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/08/13';
video_name = 'E:\06\12';
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/06/06'; % the name of the stream file (video or a collection of images)
% if the connection to the camera is live, this has to be used. 
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images'; % the name of the stream file (video or a collection of images)

% for all the time series objects and caliobration models. To use like
% this, please add Utils folder with all the subfolders to the path. Data
% subfolder keeps the timeseries objectsthank
load 'timeSeriesGHI.mat';
load 'timeSeriesPV.mat';
load 'timeSeriesTemp.mat';
load 'calib_model_01_2014_v2.mat';
load 'extrinsic_calib_model_tuned.mat';
%load 'extrinsic_calib_model_01_2014.mat';
load 'image_mask_01_2014.mat';
load 'model3D_Bigger.mat';


sunPattern = imread('sunPattern3.jpg');
modelImage = imread(model3D.fileName);

if(~exist('ocam_model'))
    display('No calibration model. Terminating...');
    return;
end

live = 0; % the pictures are live or not
track = 1;
% here we run the code 
multiObjectTrackingCV_CamCalibOcclusion(video_name,0,1,ocam_model,ecam_model,model3D,image_mask,sunPattern,modelImage,live,track,timeSeriesGHI,timeSeriesPV,timeSeriesTemp);

%%
%load 'calib_modelGood.mat';
%load 'extrinsic_calib_model.mat';