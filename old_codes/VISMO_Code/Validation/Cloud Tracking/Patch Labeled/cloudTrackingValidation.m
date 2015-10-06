%% load the required boxes
current_folder = pwd;
cd(save_folder);
load('details.mat');
load('labels.mat');
cd(current_folder);
%% load the labels before this point
for i = length(im_labelers):-1:1
    
    for j = 1:length(im_labelers(i).im_labels)
        l_str = im_labelers(i).im_labels(j);
        if(j == 1)
            mask = l_str.mask;
        else
            mask = mask|l_str.mask;
        end
        im_c = imread(l_str.name);
        bboxes(j,:) = l_str.bbox;
    end
    bboxes_all{i} = bboxes;
    bboxes = [];
    %masks{i} = repmat(mask,[1,1,3]);
end

%im = im_arr{1,1};
%im(~masks{1})=0;
%imshow(im);

%%

%centroids2 = bboxes_all{1}(:,1:2)+(bboxes_all{1}(:,3:4)/2);
%contours = cellfun(@(x) num2cell(reshape(x',2,4)',2)', num2cell([bboxes_all{1}(:,1:2),bboxes_all{1}(:,1:2)+[bboxes_all{1}(:,3),zeros(size(bboxes_all{1}(:,3)))],...
%            bboxes_all{1}(:,1:2)+bboxes_all{1}(:,3:4),bboxes_all{1}(:,1:2)+[zeros(size(bboxes_all{1}(:,4))),bboxes_all{1}(:,4)]],2)','UniformOutput',false);
%imshow(cv.drawContours(im,contours));

%%
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/05/14'; % the name of the stream file (video or a collection of images)
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/06/28';
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/08/12';
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/08/13';
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/06/06'; % the name of the stream file (video or a collection of images)
% if the connection to the camera is live, this has to be used. 
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images'; % the name of the stream file (video or a collection of images)
video_name = im_str.im_dir;

% for all the time series objects and caliobration models. To use like
% this, please add Utils folder with all the subfolders to the path. Data
% subfolder keeps the timeseries objectsthank
load 'timeSeriesGHI.mat';
load 'timeSeriesPV.mat';
load 'timeSeriesTemp.mat';
load 'calib_modelGood.mat';
load 'extrinsic_calib_model.mat';
load 'image_mask2.mat';
sunPattern = imread('sunPattern.jpg');

if(~exist('ocam_model'))
    display('No calibration model. Terminating...');
    return;
end
live = 0; % the pictures are live or not
track = 1;
% here we run the code 
errors = multiObjectTrackingCV_CamMeas(video_name,0,1,ocam_model,ecam_model,image_mask,sunPattern,live,track,timeSeriesGHI,timeSeriesPV,timeSeriesTemp,bboxes_all);

%% save the results you got
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

save('TrackingPerformance', 'errors');
save('Details', 'im_str');
save('Labels', 'im_labelers');
save('LabelsFolder','save_folder');
cd ../..
%% plotting the errors for different horizons
figure, hold on;
count = 1;
for i = 9:10:size(errors,2)
    errT = cat(1,errors(:,i).errCentroids);
    errTN = sqrt(sum(abs(errT).^2,2));
    plot(1:length(errTN),errTN,'Color',rand(1,3),'LineWidth',2);
    legendStr(count) = {strcat('Predictions for horizon of: ', num2str(i*3/60),' minutes.')};
    count = count+1;
end
title('Errors for multiple horizons (in pixels on the image plane).');
xlabel('Frame Number for the Tracking Algorithm');
ylabel('Error (in pixels on the image plane)');
legend(legendStr);

figure, hold on;
count = 1;
for i = 9:10:size(errors,2)
    errT = cat(1,errors(:,i).errCentroidsProj);
    errTN = sqrt(sum(abs(errT).^2,2));
    plot(1:length(errTN),errTN,'Color',rand(1,3),'LineWidth',2);
    legendStr(count) = {strcat('Projected predictions for horizon of: ', num2str(i*3/60),' minutes.')};
    count = count+1;
end
title('Errors for multiple horizons (in meters on the projected plane for the assumed CBH of 300 meters).');
xlabel('Frame Number for the Tracking Algorithm');
ylabel('Error (in meters on the projected plane for the assumed CBH of 300 meters)');
legend(legendStr);