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
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/06/06'; % the name of the stream file (video or a collection of images)
% if the connection to the camera is live, this has to be used. 
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images'; % the name of the stream file (video or a collection of images)
pathOld = getenv('PATH');
ctfPath = ctfroot;
opencvpath = '\opencv\build\x64\vc10\bin';
%construct and set new path
pathNew = [pathOld ';' ctfPath opencvpath];
setenv('PATH', pathNew);
% for all the time series objects and caliobration models. To use like
% this, please add Utils folder with all the subfolders to the path. Data
% subfolder keeps the timeseries objectsthank
%load 'timeSeriesGHI.mat';
%load 'timeSeriesPV.mat';
%load 'timeSeriesTemp.mat';
load 'calib_model_01_2014_v2.mat';
%load 'extrinsic_calib_model_01_2014.mat';
load 'extrinsic_calib_model_08_2013.mat';
%load 'extrinsic_calib_model_06_2014.mat';
load 'image_mask_01_2014.mat';
load 'model3D_Bigger.mat';
load 'weatherData.mat';
sunPattern = imread('sunPattern3.jpg');
modelImage = imread(model3D.fileName);

%________________________________________________________________________________________________________________
load 'PVGHIData_plus_CSK.mat';
selected_month = 8;
selected_year =  2013;
selected_day = 12;

idx=[];

for i=1:numel(GHIData)
     [y, m, d] = datevec (GHIData{i}.DateNumber);
      if selected_month == m && selected_year == y && selected_day == d      
          idx = i;      
      end     
end
%2013-11-15 10:25:00
d=datenum([2013,11,26,9,30,30])-GHIData{idx}.DateNumber;
offset = datenum([0,0,0,0,33,2])-datenum([0,0,0,0,0,d*12.3]);
%time = GHIData{idx}.Time+offset;
%occlusion_data = GHIData{idx}.GHI(3:end) - 2*GHIData{idx}.GHI(2:end-1) + GHIData{idx}.GHI(1:end-2);
%occlusion_time = GHIData{idx}.Time(3:end)+GHIData{idx}.DateNumber+offset;
timeSeriesGHI = timeseries(GHIData{idx}.GHI(:),datestr(GHIData{idx}.Time+GHIData{idx}.DateNumber+offset));
timeSeriesGHI = setinterpmethod(timeSeriesGHI,'zoh');
times_ss = pvl_maketimestruct(GHIData{idx}.Time+GHIData{idx}.DateNumber+offset, model3D.UTC);
ghi_cs = pvl_clearsky_ineichen(times_ss,model3D.SensorLocation);
timeSeriesGHI_cs = timeseries(ghi_cs,datestr(GHIData{idx}.Time+GHIData{idx}.DateNumber+offset));
timeSeriesGHI_cs = setinterpmethod(timeSeriesGHI_cs,'zoh');
%timeSeriesGHI_cs = timeseries(GHIData{idx}.ClearSkyGHI_ineichen(:),datestr(GHIData{idx}.Time+GHIData{idx}.DateNumber+offset));
%timeSeriesGHI_cs = setinterpmethod(timeSeriesGHI_cs,'zoh');
timeSeriesOcc = timeSeriesGHI;
%timeSeriesOcc = resample(timeSeriesGHI, datestr([GHIData{idx}.Time(1):datenum([0,0,0,0,0,3]):GHIData{idx}.Time(end)]+GHIData{idx}.DateNumber+offset));
occlusion_data = timeSeriesOcc.Data(3:end) - 2*timeSeriesOcc.Data(2:end-1) + timeSeriesOcc.Data(1:end-2);
occlusion_time = timeSeriesOcc.Time(3:end);
occlusion_data2 = timeSeriesOcc.Data(2:end) - timeSeriesOcc.Data(1:end-1);
occlusion_time2 = timeSeriesOcc.Time(2:end);
occlusion_data3 = occlusion_data2(3:end) - 2*occlusion_data2(2:end-1) + occlusion_data2(1:end-2);
occlusion_time3 = occlusion_time2(3:end);
%timeSeriesOcc = timeseries(occlusion_data>10 | occlusion_data2(2:end)<-10,datestr(occlusion_time+datenum(timeSeriesOcc.TimeInfo.StartDate)));
timeSeriesOcc = timeseries([occlusion_data,occlusion_data2(2:end)],datestr(occlusion_time+datenum(timeSeriesOcc.TimeInfo.StartDate)));
timeSeriesOcc = setinterpmethod(timeSeriesOcc,'zoh');
%timeSeriesOcc.DataInfo.Interpolation = tsdata.interpolation('linear');
%[ax,h1,h2]=plotyy(occlusion_time,occlusion_data,timeSeriesGHI.Time,timeSeriesGHI.Data);
%hold(ax(1),'on');
%plot(ax(1),occlusion_time2,occlusion_data2,'r');
%plotyy(occlusion_time2,occlusion_data2,timeSeriesGHI.Time,timeSeriesGHI.Data);
%plotyy(occlusion_time3,occlusion_data3,timeSeriesGHI.Time,timeSeriesGHI.Data);
%plotyy(timeSeriesOcc.Time,timeSeriesOcc.Data,timeSeriesGHI.Time,timeSeriesGHI.Data);
timeSeriesPV = timeseries(GHIData{idx}.PV_output(:),datestr(GHIData{idx}.Time+GHIData{idx}.DateNumber+offset));
timeSeriesPV = setinterpmethod(timeSeriesPV,'zoh');
timeSeriesTemp = timeseries(GHIData{idx}.ambTemp(:),datestr(GHIData{idx}.Time+GHIData{idx}.DateNumber+offset));
timeSeriesTemp = setinterpmethod(timeSeriesTemp,'zoh');
timeSeriesCBH = timeseries(((weatherData(:,2)-weatherData(:,3))/4.4)*304.8,datestr(weatherData(:,1)));
timeSeriesCBH = setinterpmethod(timeSeriesCBH,'zoh');
%_________________________________________________________________________________________________________________

if(~exist('ocam_model'))
    display('No calibration model. Terminating...');
    return;
end
%live = 0; % the pictures are live or not
%track = 0;
%record = 0;
%verbose = 1;
vid_r = 0;
vid_n = '';
% here we run the code 
fileID = fopen('conf.dat');
in_params = textscan(fileID,'%q\t%d\t%d\t%d\t%d');
video_name = in_params{1}{1};
live = in_params{2};
predict = in_params{3};
record = in_params{4};
verbose = in_params{5};
predict = 1;
verbose = 1;

%[~,~,perfs]=multiObjectTrackingCV_CamCalib('videoName',vid_n, 'recordVideo',vid_r, 'imageFolder',video_name, 'record',record, 'verbose',verbose,...
%                                           'calibModel',ocam_model, 'ecalibModel',ecam_model, '3DModel',model3D, 'mask',image_mask, 'sunPattern',sunPattern,...
%                                           'mapImage',modelImage, 'live',live, 'predict',predict, 'SeriesGHI',timeSeriesGHI, 'SeriesGHIcs',timeSeriesGHI_cs,...
%                                           'SeriesPV',timeSeriesPV, 'SeriesTemperature',timeSeriesTemp, 'SeriesOcclusion',timeSeriesOcc, 'SeriesCBH',timeSeriesCBH,...
%                                           'thresholdRGB',thRGB, 'thresholdRBR',thRBR, 'thresholdGBR',thGBR);

[~,perfs]=multiObjectTrackingCV_class(vid_n,vid_r,video_name,record,verbose,ocam_model,ecam_model,model3D,image_mask,sunPattern,modelImage,live,predict,timeSeriesGHI,timeSeriesGHI_cs,timeSeriesPV,timeSeriesTemp,timeSeriesOcc,timeSeriesCBH);
