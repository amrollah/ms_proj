%% Demo of cloud tracking code using camera inputs and calibration.
%
% Multiple object tracking with kalman filter for could tracking:
%Inputs:
% image_folder: path to the input images folder
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
close all; clear all; clc;

image_folder = 'C:\data\cav\calib';
save_name = 'extrinsic_calib_model';
% for all the time series objects and caliobration models. To use like
% this, please add Utils folder with all the subfolders to the path. Data
% subfolder keeps the timeseries objectsthank
load 'calib_model_01_2014_v2.mat';
load 'extrinsic_calib_model6.mat';
load 'image_mask_01_2014.mat';
load 'model3D_Bigger.mat';
model3D.Location.latitude = 43.5654;
model3D.Location.longitude = 11.5373;
model3D.Location.altitude = 134; % 442

sunPattern = imread('sunPattern3.jpg');

if(~exist('ocam_model'))
    display('No calibration model. Terminating...');
    return;
end
live = 0; % the pictures are live or not
frameRate = 1;
% here we run the code 
[sunPositionDetected,sunPositionReal,matchQual]=extrinsicCameraCalibration(image_folder,ocam_model,ecam_model,model3D,image_mask,sunPattern,live,frameRate);

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
cd ../

ecam_model.theta = NaN;
ecam_model.R = R;
ecam_model.Rinv = pinv(R);
save(save_name, 'ecam_model');
cd ../

%% validation
image_folder = 'C:\data\cav\calib_eval';
load 'calib_model_01_2014_v2.mat';
%load 'extrinsic_calib_model6.mat';
load 'image_mask_01_2014.mat';
%load 'extrinsic_calib_model_06_2014.mat';

cbh_given = 300;
[pix_err1,pix_err2,pix_errR1,pix_errR2] = extrinsicCameraCalibrationValidation(image_folder,ocam_model,ecam_model,model3D,image_mask,sunPattern,live,frameRate,R,cbh_given);
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