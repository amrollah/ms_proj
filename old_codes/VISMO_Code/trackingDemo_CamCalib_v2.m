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

% Input and output mechanisms of the cloud tracking code are changed, and
% the new style is adopted here. 
% The values for the required parameters are read from an .xml based
% configuration file, and the initializations are done accordingly. 

% Burak Zeydan, 22. Aug. 2014
% burak.zeydan@epfl.ch
addpath(genpath('../mexopencv-master')); %MATLAB within ABB gives error with the default path. bear this with me each time..
addpath(genpath('Utils'));

iptsetpref('ImshowBorder','tight'); %This will remove the borders at the edges of figures from MATLAB figures when opened with imshow()


str = xml2struct('config.xml');

pathOld = getenv('PATH');
ctfPath = ctfroot;
opencvpath = 'C:/opencv/build/x64/vc12/bin';
% alternateOpenCVPath = '\opencv\build\x64\vc10\bin';
%construct and set new path
if(isfield(str.properties.data,'path'))
    pathNew = [pathOld ';' ctfPath opencvpath ';' str.properties.data.path.Text];
else
    pathNew = [pathOld ';' ctfPath opencvpath];
end
setenv('PATH', pathNew);




% for all the time series objects and caliobration models. To use like
% this, please add Utils folder with all the subfolders to the path. Data
% subfolder keeps the timeseries objects
load 'calib_model_01_2014_v2.mat';
load 'extrinsic_calib_model_08_2013.mat';
load 'image_mask_01_2014.mat';
load 'model3D_Bigger.mat';
load 'weatherData.mat';
load 'PVGHIData_plus_CSK.mat';
sunPattern = imread('sunPattern3.jpg');

CHCRC_real = 'Roof-Top-Actual.jpg';
%%%%%Create a real-sized image of the map %%%% Should be done just once for
%%%%%each machine%%%% not needed anymore
% c = model3D.all; cols=-c(1)+c(2); rows=-c(3)+c(4);
% imwrite(imresize(imread(model3D.fileName),[rows cols]), CHCRC_real);
% % %%%%%%%%%


% %%%%%%% Save the Measured GHI data to text file %%%%%%
% fileName = 'C:/Users/CHFIOEZ/programs/xampp/htdocs/plotData/GHI.tsv';
% dayOfMonth = 12;
% timeStamp = datestr(GHIData{dayOfMonth}.Time,13);
% ghi = GHIData{dayOfMonth}.GHI;
% ghi2 = GHIData{dayOfMonth}.ClearSkyGHI_ineichen;
% saveTSV({timeStamp, ghi, ghi2}, fileName);
%%%%%%%%%%%

if(isfield(str.properties.data,'calibModel'))
    load(str.properties.data.calibModel.Text);
end
if(isfield(str.properties.data,'ecalibModel'))
    load(str.properties.data.ecalibModel.Text);
end
if(isfield(str.properties.data,'imageMask'))
    load(str.properties.data.imageMask.Text);
end
if(isfield(str.properties.data,'Model3D'))
    load(str.properties.data.Model3D.Text);
end
if(isfield(str.properties.data,'weatherData'))
    load(str.properties.data.weatherData.Text);
end
if(isfield(str.properties.data,'GHIData'))
    load(str.properties.data.GHIData.Text);
end
if(isfield(str.properties.data,'sunPattern'))
    sunPattern = imread(str.properties.data.sunPattern.Text);
end

modelImage = imread(model3D.fileName);
model3D.fileNameActual = CHCRC_real;
%________________________________________________________________________________________________________________
tic
timesM = [];
dataGHI = [];
dataPV = [];
dataTemp = [];
tic
for i=1:numel(GHIData)
    d=datenum([2013,11,26,9,30,30])-GHIData{i}.DateNumber;
    offset = datenum([0,0,0,0,33,2])-datenum([0,0,0,0,0,d*12.3]);
    timesM = [timesM;GHIData{i}.Time+GHIData{i}.DateNumber+offset];
    dataGHI = [dataGHI;GHIData{i}.GHI(:)];
    dataPV = [dataPV;GHIData{i}.PV_output(:)];
    dataTemp = [dataTemp;GHIData{i}.ambTemp(:)];
end
toc
timesMStr = datestr(timesM);
timeSeriesGHI = timeseries(dataGHI,timesMStr);
timeSeriesGHI = setinterpmethod(timeSeriesGHI,'zoh');
timeSeriesPV = timeseries(dataGHI,timesMStr);
timeSeriesPV = setinterpmethod(timeSeriesPV,'zoh');
timeSeriesTemp = timeseries(dataGHI,timesMStr);
timeSeriesTemp = setinterpmethod(timeSeriesTemp,'zoh');

times_ss = pvl_maketimestruct(timesM, model3D.UTC);
ghi_cs = pvl_clearsky_ineichen(times_ss,model3D.SensorLocation);
timeSeriesGHI_cs = timeseries(ghi_cs,timesMStr);
timeSeriesGHI_cs = setinterpmethod(timeSeriesGHI_cs,'zoh');

timeSeriesOcc = timeSeriesGHI;
occlusion_data = timeSeriesOcc.Data(3:end) - 2*timeSeriesOcc.Data(2:end-1) + timeSeriesOcc.Data(1:end-2);
occlusion_time = timeSeriesOcc.Time(3:end);
occlusion_data2 = timeSeriesOcc.Data(2:end) - timeSeriesOcc.Data(1:end-1);
occlusion_time2 = timeSeriesOcc.Time(2:end);
occlusion_data3 = occlusion_data2(3:end) - 2*occlusion_data2(2:end-1) + occlusion_data2(1:end-2);
occlusion_time3 = occlusion_time2(3:end);
timeSeriesOcc = timeseries([occlusion_data,occlusion_data2(2:end)],datestr(occlusion_time+datenum(timeSeriesOcc.TimeInfo.StartDate)));
timeSeriesOcc = setinterpmethod(timeSeriesOcc,'zoh');

timeSeriesCBH = timeseries(((weatherData(:,2)-weatherData(:,3))/4.4)*304.8,datestr(weatherData(:,1)));
timeSeriesCBH = setinterpmethod(timeSeriesCBH,'zoh');
toc
%_________________________________________________________________________________________________________________

if(~exist('ocam_model') || ~exist('ecam_model'))
    display('Missing calibration model. Terminating...');
    return;
end

% Commenting out Burak's dir.
% imageFolder = 'C:/Users/chbuzey/Documents/VisualStudio2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/08/13';
% imageFolder = 'C:/Users/CHFIOEZ/Documents/git/cloudtracking/codes/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/08/13'; %note that currently the code only takes config.xml into consideration.


% if the connection to the camera is live, the below one should be used. 
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images'; % the name of the stream file (video or a collection of images)
%video_name = ' C:/Users/CHFIOEZ/Documents/git/cloudtracking/codes/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images'; %TODO: check this

live = 0; % the pictures are live or not
predict = 1;
record = 0;
verbose = 1;
recordVideo = 0;
videoName = '';
horizon = 25;
frameRate = 4; %in terms of seconds per frame.
rgbThreshold = 100;
rbrThreshold = 0.72;
gbrThreshold = 0.8;
circumsolarSize = 200;
suneliminationSize = 80;

if(verbose)
    fprintf('Camera frame rate is %ffps; one frame captured per %dseconds \n', 1/(3*frameRate), 3*frameRate);
end

if(isfield(str.properties.algorithm,'live'))
    live = str2num(str.properties.algorithm.live.Text);
end
if(isfield(str.properties.algorithm,'predict'))
    predict = str2num(str.properties.algorithm.predict.Text);
end
if(isfield(str.properties.algorithm,'record'))
    record = str2num(str.properties.algorithm.record.Text);
end
if(isfield(str.properties.algorithm,'verbose'))
    verbose = str2num(str.properties.algorithm.verbose.Text);
end
if(isfield(str.properties.algorithm,'recordVideo'))
    recordVideo = str2num(str.properties.algorithm.recordVideo.Text);
end
if(isfield(str.properties.algorithm,'videoName'))
    videoName = str.properties.algorithm.videoName.Text;
end
if(isfield(str.properties.algorithm,'imageFolder'))
    imageFolder = str.properties.algorithm.imageFolder.Text;
end
if(isfield(str.properties.algorithm,'saveImages'))
    saveImages = str2num(str.properties.algorithm.saveImages.Text);
end
if(isfield(str.properties.algorithm,'pathRelativeSaveImages'))
    pathRelativeSaveImages = str.properties.algorithm.pathRelativeSaveImages.Text;
end



if(isfield(str.properties.parameters,'horizon'))
    horizon = str2num(str.properties.parameters.horizon.Text);
end
if(isfield(str.properties.parameters,'frameRate'))
    frameRate = str2num(str.properties.parameters.frameRate.Text);
end
if(isfield(str.properties.parameters,'rgbThreshold'))
    rgbThreshold = str2num(str.properties.parameters.rgbThreshold.Text);
end
if(isfield(str.properties.parameters,'rbrThreshold'))
    rbrThreshold = str2num(str.properties.parameters.rbrThreshold.Text);
end
if(isfield(str.properties.parameters,'gbrThreshold'))
    gbrThreshold = str2num(str.properties.parameters.gbrThreshold.Text);
end
if(isfield(str.properties.parameters,'circumsolarSize'))
    circumsolarSize = str2num(str.properties.parameters.circumsolarSize.Text);
end
if(isfield(str.properties.parameters,'suneliminationSize'))
    suneliminationSize = str2num(str.properties.parameters.suneliminationSize.Text);
end

[~,~,perfs]=multiObjectTrackingCV_CamCalib('videoName',videoName, 'recordVideo',recordVideo, 'imageFolder',imageFolder, 'record',record, 'verbose',verbose,...
                                           'calibModel',ocam_model, 'ecalibModel',ecam_model, '3DModel',model3D, 'mask',image_mask, 'sunPattern',sunPattern,...
                                           'mapImage',modelImage, 'live',live, 'predict',predict, 'SeriesGHI',timeSeriesGHI, 'SeriesGHIcs',timeSeriesGHI_cs,...
                                           'SeriesPV',timeSeriesPV, 'SeriesTemperature',timeSeriesTemp, 'SeriesOcclusion',timeSeriesOcc, 'SeriesCBH',timeSeriesCBH,...
                                           'horizon',horizon, 'frameRate',frameRate, 'thresholdRGB', rgbThreshold, 'thresholdRBR',rbrThreshold, 'thresholdGBR',gbrThreshold,...
                                           'circumsolarSize',circumsolarSize, 'suneliminationSize',suneliminationSize, 'saveImages', saveImages, ...
                                           'pathRelativeSaveImages', pathRelativeSaveImages);

%[~,perfs]=multiObjectTrackingCV_class(vid_n,vid_r,video_name,record,verbose,ocam_model,ecam_model,model3D,image_mask,sunPattern,modelImage,live,track,timeSeriesGHI,timeSeriesGHI_cs,timeSeriesPV,timeSeriesTemp,timeSeriesOcc,timeSeriesCBH);
