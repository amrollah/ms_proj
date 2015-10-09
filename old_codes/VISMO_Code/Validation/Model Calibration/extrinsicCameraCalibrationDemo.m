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
addpath(genpath('../../../mexopencv-master')); %MATLAB within ABB gives error with the default path. bear this with me each time..
addpath(genpath('../../Utils'));

%video_name = 'C:\Users\chbuzey\Documents\MATLAB\Archive\Validation Sequences\Model Validation\1';
%video_name = 'E:\06\calib12';
%video_name = 'E:\VISMO IMAGES\2013\08\ecalib_set';
video_name = 'C:\data\04';
save_name = 'extrinsic_calib_model_08_2013';
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/07/21';
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/06/06'; % the name of the stream file (video or a collection of images)
% if the connection to the camera is live, this has to be used. 
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images'; % the name of the stream file (video or a collection of images)

% for all the time series objects and caliobration models. To use like
% this, please add Utils folder with all the subfolders to the path. Data
% subfolder keeps the timeseries objectsthank
load 'calib_model_01_2014_v2.mat';
load 'extrinsic_calib_model6.mat';
load 'image_mask_01_2014.mat';
load 'model3D_Bigger.mat';
%load 'image_mask2.mat';

sunPattern = imread('sunPattern3.jpg');

if(~exist('ocam_model'))
    display('No calibration model. Terminating...');
    return;
end
live = 0; % the pictures are live or not
frameRate = 20;
% here we run the code 
[sunPositionDetected,sunPositionReal,matchQual]=extrinsicCameraCalibration(video_name,ocam_model,ecam_model,model3D,image_mask,sunPattern,live,frameRate);

%% optimization
X_c = sunPositionDetected';
X_w = sunPositionReal';
X_w(:,3) = -X_w(:,3);
R = zeros(3,3);
%X_c = X_c(10:end-10,:);
%X_w = X_w(10:end-10,:);
%X_c = [X_c,ones(size(X_c,1),1)];
X = blkdiag(X_c,X_c,X_c);
y = X_w(:);
B = inv(X'*X)*X'*y;
R = reshape(B,3,length(B)/3)';

%for i = 1:3
%    y = X_w(:,i);
%    B = inv(X_c'*X_c)*X_c'*y;
%    R(i,:) = B';
%end

%R = reshape(fmincon(@(x)(norm(blkdiag(X_c,X_c,X_c)*x-X_w(:))^2),ones(9,1),[],[],[],[],[],[],@(x)deal(0,x'*x-3)),3,3)';
%R = fmincon(@(x)(norm(x*X_c'-X_w')^2),eye(3,3),[],[],[],[],[],[],@(x)deal([],x'*x-eye(3,3)));
figure
scatter3(X_w(:,1),X_w(:,2),-X_w(:,3),9,'sb','LineWidth',4);
hold on;
scatter3(X_c(:,1),X_c(:,2),-X_c(:,3),9,'or','LineWidth',4);
%plot3(X_w(:,1),X_w(:,2),-X_w(:,3),'g');
%plot3(X_c(:,1),X_c(:,2),-X_c(:,3),'k');
XX = (R*X_c')';
scatter3(XX(:,1),XX(:,2),-XX(:,3),9,'xy','LineWidth',4);
%% save values
folderName = 'Results';

e_val = exist(folderName);
if ~(e_val == 7)
    mkdir(folderName);
end
cd(folderName);

c_list = ls;
b_num = sort(str2num(c_list(3:end,:)));

if isempty(b_num)
    b_num = 0;
else
    b_num = b_num(end);
end

b_num = b_num + 1;
b_str = num2str(b_num);
mkdir(b_str);
cd(b_str);

save('Rotation_Matrix','R');
save('SunDetections','sunPositionDetected');
save('SunReal','sunPositionReal');
cd ../..

ecam_model.theta = NaN;
ecam_model.R = R;
ecam_model.Rinv = pinv(R);
cd ../../Utils/Data
save(save_name, 'ecam_model');
cd '../../Validation/Model Calibration'
%% test results
%sunPattern = imread('sunPattern3.jpg');
%video_name = 'C:\Users\chbuzey\Documents\MATLAB\Archive\Validation Sequences\Model Validation\1';
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/08/13';
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/08/12';
%video_name = 'E:\VISMO IMAGES\2013\07\29';
%video_name = 'F:\VISMO_img\2013\07\23';
%frameRate = 50;
%extrinsicCameraCalibrationValidation(video_name,ocam_model,ecam_model,model3D,image_mask,sunPattern,live,frameRate,R);
%video_name = 'C:\Users\chbuzey\Documents\MATLAB\Archive\Validation Sequences\Model Validation\1';
%video_name = 'E:\06\calib12';
video_name = 'C:\data\04_validation';
load 'calib_model_01_2014_v2.mat';
load 'extrinsic_calib_model_06_2014.mat';
%load 'extrinsic_calib_model_tuned.mat';
load 'image_mask_01_2014.mat';
R = ecam_model.R;
load 'extrinsic_calib_model6.mat';
load 'model3D_Bigger.mat';
sunPattern = imread('sunPattern3.jpg');
live = 0;
frameRate = 50;
cbh_given = 1500;
[pix_err1,pix_err2,pix_errR1,pix_errR2] = extrinsicCameraCalibrationValidation(video_name,ocam_model,ecam_model,model3D,image_mask,sunPattern,live,frameRate,R,cbh_given);
disp(['Previous Backprojection Error: ',num2str(mean(pix_err2)),' ',setstr(177),' ', num2str(std(pix_err2)),' pixels.']);
disp(['Backprojection Error: ',num2str(mean(pix_err1)),' ',setstr(177),' ', num2str(std(pix_err1)),' pixels.']);
disp(['Previous Projection Error at a 1500 meter height: ',num2str(mean(pix_errR2)),' ',setstr(177),' ', num2str(std(pix_errR2)),' meters.']);
disp(['Projection Error at a 1500 meter height: ',num2str(mean(pix_errR1)),' ',setstr(177),' ', num2str(std(pix_errR1)),' meters.']);
figure, plot([1:size(pix_err1,2)],pix_err1,[1:size(pix_err2,2)],pix_err2);
figure, plot([1:size(pix_err1,2)],pix_err1);
figure, plot([1:size(pix_err2,2)],pix_err2);
figure, plot([1:size(pix_errR1,2)],pix_errR1,[1:size(pix_errR2,2)],pix_errR2);
figure, plot([1:size(pix_errR1,2)],pix_errR1);
figure, plot([1:size(pix_errR2,2)],pix_errR2);

%%
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/08/12';

load 'calib_model_01_2014_v2.mat';
load 'extrinsic_calib_model_01_2014.mat';
%load 'extrinsic_calib_model_06_2014.mat';
%load 'extrinsic_calib_model_tuned.mat';
load 'image_mask_01_2014.mat';
R = ecam_model.R;
load 'extrinsic_calib_model6.mat';
load 'model3D_Bigger.mat';
sunPattern = imread('sunPattern3.jpg');

deg2radf = pi/180;
DN = datenum([2013,8,12,7,0,1]);
Time = pvl_maketimestruct(DN, model3D.UTC);
[SunAz, ~, ApparentSunEl, ~]=pvl_ephemeris(Time, model3D.Location);
zenith = 90-ApparentSunEl;
azimuth = SunAz;
[x,y,z] = sph2cart((90-azimuth)*deg2radf,(90-zenith)*deg2radf,1);
sunPositionReal = [x;y;z];
pinv(R)*[x;y;-z]

image_mask2 = ones(ocam_model.height,ocam_model.width);
[X,Y] = meshgrid(1:ocam_model.height,1:ocam_model.width);
P = [X(:)';Y(:)'];
p = cam2world(P,ocam_model);
p = R*p;
maskTemp = p(3,:)<-0.05;
indicesTemp = sub2ind([ocam_model.height,ocam_model.width],P(1,:),P(2,:));
image_mask2(indicesTemp) = maskTemp;
image_mask2 = logical(repmat(image_mask2,[1,1,3]));
image_mask = image_mask2;

frameRate = 50;
live = 0;
cbh_given = 1500;
[pix_err1,pix_err2,pix_errR1,pix_errR2] = extrinsicCameraCalibrationAdaptive(video_name,ocam_model,ecam_model,model3D,image_mask,sunPattern,live,frameRate,R,cbh_given);
disp(['Previous Backprojection Error: ',num2str(mean(pix_err2)),' ',setstr(177),' ', num2str(std(pix_err2)),' pixels.']);
disp(['Backprojection Error: ',num2str(mean(pix_err1)),' ',setstr(177),' ', num2str(std(pix_err1)),' pixels.']);
disp(['Previous Projection Error at a 1500 meter height: ',num2str(mean(pix_errR2)),' ',setstr(177),' ', num2str(std(pix_errR2)),' meters.']);
disp(['Projection Error at a 1500 meter height: ',num2str(mean(pix_errR1)),' ',setstr(177),' ', num2str(std(pix_errR1)),' meters.']);
figure, plot([1:size(pix_err1,2)],pix_err1,[1:size(pix_err2,2)],pix_err2);
figure, plot([1:size(pix_err1,2)],pix_err1);
figure, plot([1:size(pix_err2,2)],pix_err2);
figure, plot([1:size(pix_errR1,2)],pix_errR1,[1:size(pix_errR2,2)],pix_errR2);
figure, plot([1:size(pix_errR1,2)],pix_errR1);
figure, plot([1:size(pix_errR2,2)],pix_errR2);


%% for changing CBH how does the external calibration error change...
load 'calib_model_01_2014_v2.mat';
%load 'extrinsic_calib_model_01_2014.mat';
load 'extrinsic_calib_model_tuned.mat';
load 'image_mask_01_2014.mat';
R = ecam_model.R;
load 'extrinsic_calib_model6.mat';
load 'model3D_Bigger.mat';
sunPattern = imread('sunPattern3.jpg');
frameRate = 50;
live = 0;

cbh_base = 50;
cbh_low = 1;
cbh_up = 5;
cbh_scale = 100/cbh_up;
for i = cbh_low:cbh_up
    cbh_given = cbh_base * i * cbh_scale;
    video_name = 'C:\data\12';
    [pix_err1,pix_err2,pix_errR1,pix_errR2,pix_errR3,pix_errR4,pix_errR5,pix_errR6] = extrinsicCameraCalibrationValidation(video_name,ocam_model,ecam_model,model3D,image_mask,sunPattern,live,frameRate,R,cbh_given);
    mean_bb_err_pix1(i) = mean(pix_err1);
    mean_bb_err_pix2(i) = mean(pix_err2);
    max_bb_err_pix1(i) = max(pix_err1);
    max_bb_err_pix2(i) = max(pix_err2);
    min_bb_err_pix1(i) = min(pix_err1);
    min_bb_err_pix2(i) = min(pix_err2);
    med_bb_err_pix1(i) = median(pix_err1);
    med_bb_err_pix2(i) = median(pix_err2);
    mean_bb_err_m1(i) = mean(pix_errR1);
    mean_bb_err_m2(i) = mean(pix_errR2);
    max_bb_err_m1(i) = max(pix_errR1);
    max_bb_err_m2(i) = max(pix_errR2);
    min_bb_err_m1(i) = min(pix_errR1);
    min_bb_err_m2(i) = min(pix_errR2);
    med_bb_err_m1(i) = median(pix_errR1);
    med_bb_err_m2(i) = median(pix_errR2);
    mean_bb_err_m3(i) = mean(pix_errR3);
    mean_bb_err_m4(i) = mean(pix_errR4);
    max_bb_err_m3(i) = max(pix_errR3);
    max_bb_err_m4(i) = max(pix_errR4);
    min_bb_err_m3(i) = min(pix_errR3);
    min_bb_err_m4(i) = min(pix_errR4);
    med_bb_err_m3(i) = median(pix_errR3);
     med_bb_err_m4(i) = median(pix_errR4);
    mean_bb_err_m5(i) = mean(pix_errR5);
    mean_bb_err_m6(i) = mean(pix_errR6);
    max_bb_err_m5(i) = max(pix_errR5);
    max_bb_err_m6(i) = max(pix_errR6);
    min_bb_err_m5(i) = min(pix_errR5);
    min_bb_err_m6(i) = min(pix_errR6);
    med_bb_err_m5(i) = median(pix_errR5);
    med_bb_err_m6(i) = median(pix_errR6);
    
%     video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/08/12';
%     [pix_err1,pix_err2,pix_errR1,pix_errR2] = extrinsicCameraCalibrationAdaptive(video_name,ocam_model,ecam_model,model3D,image_mask,sunPattern,live,frameRate,R,cbh_given);
%     mean_bb_err_pix3(i) = mean(pix_err1);
%     mean_bb_err_pix4(i) = mean(pix_err2);
%     max_bb_err_pix3(i) = max(pix_err1);
%     max_bb_err_pix4(i) = max(pix_err2);
%     min_bb_err_pix3(i) = min(pix_err1);
%     min_bb_err_pix4(i) = min(pix_err2);
%     med_bb_err_pix3(i) = median(pix_err1);
%     med_bb_err_pix4(i) = median(pix_err2);
%     mean_bb_err_m3(i) = mean(pix_errR1);
%     mean_bb_err_m4(i) = mean(pix_errR2);
%     max_bb_err_m3(i) = max(pix_errR1);
%     max_bb_err_m4(i) = max(pix_errR2);
%     min_bb_err_m3(i) = min(pix_errR1);
%     min_bb_err_m4(i) = min(pix_errR2);
%     med_bb_err_m3(i) = median(pix_errR1);
%     med_bb_err_m4(i) = median(pix_errR2);
%     
%     video_name = 'C:\Users\chbuzey\Documents\MATLAB\Archive\Validation Sequences\Model Validation\1';
%     [pix_err1,pix_err2,pix_errR1,pix_errR2] = extrinsicCameraCalibrationAdaptive(video_name,ocam_model,ecam_model,model3D,image_mask,sunPattern,live,frameRate,eye(3),cbh_given);
%     mean_bb_err_pix5(i) = mean(pix_err1);
%     mean_bb_err_pix6(i) = mean(pix_err2);
%     max_bb_err_pix5(i) = max(pix_err1);
%     max_bb_err_pix6(i) = max(pix_err2);
%     min_bb_err_pix5(i) = min(pix_err1);
%     min_bb_err_pix6(i) = min(pix_err2);
%     med_bb_err_pix5(i) = median(pix_err1);
%     med_bb_err_pix6(i) = median(pix_err2);
%     mean_bb_err_m5(i) = mean(pix_errR1);
%     mean_bb_err_m6(i) = mean(pix_errR2);
%     max_bb_err_m5(i) = max(pix_errR1);
%     max_bb_err_m6(i) = max(pix_errR2);
%     min_bb_err_m5(i) = min(pix_errR1);
%     min_bb_err_m6(i) = min(pix_errR2);
%     med_bb_err_m5(i) = median(pix_errR1);
%     med_bb_err_m6(i) = median(pix_errR2);
end

figure, plot([cbh_low:cbh_up]*cbh_scale*cbh_base,max_bb_err_m1,'-',[cbh_low:cbh_up]*cbh_scale*cbh_base,mean_bb_err_m1,'-',[cbh_low:cbh_up]*cbh_scale*cbh_base,med_bb_err_m1,'-',[cbh_low:cbh_up]*cbh_scale*cbh_base,min_bb_err_m1,'LineWidth',3);
legend('Maximum','Mean','Median','Minimum','Location','Northwest');
grid on;
set(gca,'FontSize',17);
xlabel('Cloud Base Height (m)','FontSize', 18);
ylabel('Projection Error (m)','FontSize', 18);
title('Projection Error vs. Cloud Base Height','FontSize', 18);

figure, plot([cbh_low:cbh_up]*cbh_scale*cbh_base,mean_bb_err_m1,[cbh_low:cbh_up]*cbh_scale*cbh_base,med_bb_err_m1,[cbh_low:cbh_up]*cbh_scale*cbh_base,min_bb_err_m1,'LineWidth',3);
legend('Mean','Median','Minimum','Location','Northwest');
grid on;
set(gca,'FontSize',17);
xlabel('Cloud Base Height (m)','FontSize', 18);
ylabel('Projection Error (m)','FontSize', 18);
title('System Error vs. Cloud Base Height','FontSize', 18);


figure, plot([cbh_low:cbh_up]*cbh_scale*cbh_base,mean_bb_err_m1+mean_err(20:20:100),[cbh_low:cbh_up]*cbh_scale*cbh_base,med_bb_err_m1+med_err(20:20:100),[cbh_low:cbh_up]*cbh_scale*cbh_base,min_bb_err_m1+min_err(20:20:100),'LineWidth',3);
legend('Mean','Median','Minimum','Location','Northwest');
grid on;
set(gca,'FontSize',17);
xlabel('Cloud Base Height (m)','FontSize', 18);
ylabel('Projection Error (m)','FontSize', 18);
title('System Error vs. Cloud Base Height','FontSize', 18);

figure, plot([cbh_low:cbh_up]*cbh_scale*cbh_base,mean_bb_err_m3+mean_err(20:20:100),[cbh_low:cbh_up]*cbh_scale*cbh_base,med_bb_err_m3+med_err(20:20:100),[cbh_low:cbh_up]*cbh_scale*cbh_base,min_bb_err_m3+min_err(20:20:100),'LineWidth',3);
legend('Mean','Median','Minimum','Location','Northwest');
grid on;
set(gca,'FontSize',17);
xlabel('Cloud Base Height (m)','FontSize', 18);
ylabel('Projection Error (m)','FontSize', 18);
title('System Error vs. Cloud Base Height','FontSize', 18);

figure, plot([cbh_low:cbh_up]*cbh_scale*cbh_base,mean_bb_err_m1,[cbh_low:cbh_up]*cbh_scale*cbh_base,med_bb_err_m1,[cbh_low:cbh_up]*cbh_scale*cbh_base,min_bb_err_m1,...
    [cbh_low:cbh_up]*cbh_scale*cbh_base,mean_bb_err_m3,[cbh_low:cbh_up]*cbh_scale*cbh_base,med_bb_err_m3,[cbh_low:cbh_up]*cbh_scale*cbh_base,min_bb_err_m3,...
    [cbh_low:cbh_up]*cbh_scale*cbh_base,mean_bb_err_m5,[cbh_low:cbh_up]*cbh_scale*cbh_base,med_bb_err_m5,[cbh_low:cbh_up]*cbh_scale*cbh_base,min_bb_err_m5,'LineWidth',3);
legend('Mean','Median','Minimum','Mean Upper','Median Upper','Minimum Upper','Mean Lower','Median Lower','Minimum Lower','Location','Northwest');
grid on;
set(gca,'FontSize',17);
xlabel('Cloud Base Height (m)','FontSize', 18);
ylabel('Projection Error (m)','FontSize', 18);
title('Projection Error vs. Cloud Base Height','FontSize', 18);

figure, plot([cbh_low:cbh_up]*cbh_scale*cbh_base,mean_bb_err_m1+mean_err(20:20:100),'-r',[cbh_low:cbh_up]*cbh_scale*cbh_base,mean_bb_err_m3+mean_err(20:20:100),'--r',...
    [cbh_low:cbh_up]*cbh_scale*cbh_base,med_bb_err_m1+med_err(20:20:100),'-b',[cbh_low:cbh_up]*cbh_scale*cbh_base,med_bb_err_m3+med_err(20:20:100),'--b',...
    [cbh_low:cbh_up]*cbh_scale*cbh_base,min_bb_err_m1+min_err(20:20:100),'-g',[cbh_low:cbh_up]*cbh_scale*cbh_base,min_bb_err_m3+min_err(20:20:100),'--g','LineWidth',3);
legend('Mean (lower)','Mean (upper)','Median (lower)','Median (upper)','Minimum (lower)','Minimum (upper)','Location','Northwest');
grid on;
set(gca,'FontSize',17);
xlabel('Cloud Base Height (m)','FontSize', 18);
ylabel('Projection Error (m)','FontSize', 18);
title('Projection Error vs. Cloud Base Height','FontSize', 18);

figure, plot([cbh_low:cbh_up]*cbh_scale*cbh_base,mean_bb_err_m1+mean_err(20:20:100),[cbh_low:cbh_up]*cbh_scale*cbh_base,med_bb_err_m1+med_err(20:20:100),[cbh_low:cbh_up]*cbh_scale*cbh_base,min_bb_err_m1+min_err(20:20:100),'LineWidth',3);
legend('Mean','Median','Minimum','Location','Northwest');
grid on;
set(gca,'FontSize',17);
xlabel('Cloud Base Height (m)','FontSize', 18);
ylabel('Projection Error (m)','FontSize', 18);
title('System Error vs. Cloud Base Height','FontSize', 18);

figure, plot([cbh_low:cbh_up]*cbh_scale*cbh_base,mean_bb_err_m5,[cbh_low:cbh_up]*cbh_scale*cbh_base,med_bb_err_m5,[cbh_low:cbh_up]*cbh_scale*cbh_base,min_bb_err_m5,'LineWidth',3);
legend('Mean','Median','Minimum','Location','Northwest');
grid on;
set(gca,'FontSize',17);
xlabel('Cloud Base Height (m)','FontSize', 18);
ylabel('Projection Error (m)','FontSize', 18);
title('Projection Error vs. Cloud Base Height','FontSize', 18);
%%
load 'calib_model_01_2014_v2.mat';
%load 'extrinsic_calib_model_01_2014.mat';
load 'extrinsic_calib_model_tuned.mat';
load 'image_mask_01_2014.mat';
R = ecam_model.R;
load 'extrinsic_calib_model6.mat';
load 'model3D_Bigger.mat';
sunPattern = imread('sunPattern3.jpg');
frameRate = 50;
live = 0;
cbh_given = 100;
video_name = 'C:\data\12';
frame_res = extrinsicCameraCalibrationAvg(video_name,ocam_model,ecam_model,model3D,image_mask,sunPattern,live,frameRate,R,cbh_given);


