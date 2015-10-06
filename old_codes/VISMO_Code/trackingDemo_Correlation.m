% B. Zeydan, 08. Oct. 2013



%video_name = 'Cloud_Tracking_Videos/slow_1/sky_low.mp4'; % the name of the stream file (video or a collection of images)
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/08/12';
%video_name = 'C:\Users\chbuzey\Documents\MATLAB\CH-CR roof images\Utils\Images\august 12';
%video_name = 'C:\Users\chbuzey\Documents\MATLAB\CH-CR roof images\Utils\Images\september 13';
video_name = 'E:\06\12';
%video_name = 'C:\Users\chbuzey\Documents\MATLAB\CH-CR roof images\Cloud_Tracking_Videos\Images\Simple_Correlation';
sunPattern = imread('sunPattern.jpg');
load 'calib_model_01_2014_v2.mat';
load 'extrinsic_calib_model_01_2014.mat';
load 'image_mask_01_2014.mat';
load 'model3D_Bigger.mat';
load 'weatherData.mat';
sunPattern = imread('sunPattern3.jpg');
modelImage = imread(model3D.fileName);
multiObjectTrackingCV_Correlation(video_name,0,1,2000,image_mask,sunPattern,ocam_model,ecam_model);