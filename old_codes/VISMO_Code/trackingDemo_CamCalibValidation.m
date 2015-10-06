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
video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/08/13';
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/06/06'; % the name of the stream file (video or a collection of images)
% if the connection to the camera is live, this has to be used.
%video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images'; % the name of the stream file (video or a collection of images)

% for all the time series objects and caliobration models. To use like
% this, please add Utils folder with all the subfolders to the path. Data
% subfolder keeps the timeseries objectsthank
load 'timeSeriesGHI.mat';
load 'timeSeriesPV.mat';
load 'timeSeriesTemp.mat';
load 'calib_model_01_2014.mat';
load 'extrinsic_calib_model_01_2014.mat';
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
record = 1;
verbose = 0;
% here we run the code

[obj_r,recorded_tracks] = multiObjectTrackingCV_CamCalib(video_name,record,verbose,ocam_model,ecam_model,model3D,image_mask,sunPattern,modelImage,live,track,timeSeriesGHI,timeSeriesPV,timeSeriesTemp);
close all
%%
errs = multiObjectTrackingCV_CamCalibValidation(obj_r,recorded_tracks);

if(record)
    e_val = exist('recordings');
    if ~(e_val == 7)
        mkdir recordings;
    end
    cd recordings;
    c_list = sort(cell2mat(cellfun(@(x) str2num(x),cellstr(ls),'UniformOutput',false)));
    
    if ~isempty(c_list)
        b_num = c_list(end);
    else
        b_num = 0;
    end
    
    b_num = b_num + 1;
    b_str = num2str(b_num);
    mkdir(b_str);
    cd(b_str);
    save('recorded_tracks', 'recorded_tracks','-v7.3');
    save('object', 'obj_r');
    save('errors', 'errs','-v7.3');
    cd ../..
    break;
end

%%
if(~record)
    max = 0;
    ind = 0;
    for i=1:size(errs,2)
        errT = errs(:,i,:);
        errTH = errT(:,:,1);
        activeMeasObs = sum(~cellfun(@isempty,{errTH.errCentroidN}));
        size_mat(i) = activeMeasObs;
        centroidObs = cat(1,errTH(:).centroidObs);
        if(activeMeasObs>1)
            ke = 0; kb = 0;
            while((sum(centroidObs(end-ke,:)==-1)==2)&&(ke<size(centroidObs,1)-1))
                ke = ke+1;
            end
            while((sum(centroidObs(1+kb,:)==-1)==2)&&(kb<size(centroidObs,1)-1))
                kb = kb+1;
            end
            if((sum(centroidObs(end-ke,:)==-1)==2) || (sum(centroidObs(1+kb,:)==-1)==2) || (1+kb > size(centroidObs,1)-ke))
                disp_mat(i) = 0;
            else
                disp_mat(i) = norm(centroidObs(1+kb,:)-centroidObs(end-ke,:));
            end
        else
            disp_mat(i) = 0;
        end
        %if(activeMeasObs>max)
        %    ind = i;
        %    max = activeMeasObs;
        %end
    end
    [ss,ssi] = sort(size_mat);
    [sd,sdi] = sort(disp_mat);
    [ind_int,ind_s,ind_d] = intersect(ssi(end-100:end),sdi(end-100:end));
    [s_int,s_ind]=sort(ind_s+ind_d);
    ind = ind_int(s_ind(end-2));
    ind = 231;
    %ind = 272;
    %ind = 266;
    %ind = 224;
    %ind = 94;
    %ind_s = find(size_mat==ss(end-20));
    %ind_d = find(disp_mat==sd(end));
    %ind_s  = ssi(end-20);
    %ind_d  = sdi(end);
    %ind = ind_d;
    errT = errs(:,ind(end),:);
    horizon = 2;
    errTH = errT(:,:,horizon);
    xdim = 1:length(errTH);
    
    centroidObs = cat(1,errTH(:).centroidObs);
    centroidPred = cat(1,errTH(:).centroidPred);
    
    activeMeasObs = ~cellfun(@isempty,{errTH.centroidObs});
    activeMeasPred = ~cellfun(@isempty,{errTH.centroidPred});
    plot(xdim(activeMeasObs),centroidObs(:,2),':gs',xdim(activeMeasPred),centroidPred(:,2),'--go','Color',rand(1,3),'LineWidth',2);
    
    %activeMeas = ~cellfun(@isempty,{errTH.centroidObs});
    %plot(xdim(activeMeas),centroidObs(:,2),,'Color',rand(1,3),'LineWidth',2);
    
    %activeMeas = ~cellfun(@isempty,{errTH.centroidPred});
    %plot(xdim(activeMeas),centroidPred(:,2),'Color',rand(1,3),'LineWidth',2);
    activeMeas = ~cellfun(@isempty,{errTH.errCentroidN});
    plot(xdim(activeMeas),[errTH(:).errCentroidN],'-x','Color',rand(1,3),'LineWidth',3,'MarkerSize',10);
    
    centroidProjObs = cat(1,errTH(:).centroidProjObs);
    centroidProjPred = cat(1,errTH(:).centroidProjPred);
    
    activeMeasObs = ~cellfun(@isempty,{errTH.centroidProjObs});
    activeMeasPred = ~cellfun(@isempty,{errTH.centroidProjPred});
    plot(xdim(activeMeasObs),centroidProjObs(:,2),':gs',xdim(activeMeasPred),centroidProjPred(:,2),'--go','Color',rand(1,3),'LineWidth',2);
    
    %activeMeas = ~cellfun(@isempty,{errTH.centroidObs});
    %plot(xdim(activeMeas),centroidObs(:,2),,'Color',rand(1,3),'LineWidth',2);
    
    %activeMeas = ~cellfun(@isempty,{errTH.centroidPred});
    %plot(xdim(activeMeas),centroidPred(:,2),'Color',rand(1,3),'LineWidth',2);
    activeMeas = ~cellfun(@isempty,{errTH.errCentroidN});
    centroidObs = cat(1,errTH(activeMeas).centroidObs);
    plot(centroidObs(:,2),centroidObs(:,1),':gs',centroidPred(:,2),centroidPred(:,1),'--go','Color',rand(1,3),'LineWidth',2);
    
    figure; hold on;
    plot(centroidObs(:,2),centroidObs(:,1),':g','Color',rand(1,3),'LineWidth',2);
    plot(centroidPred(:,2),centroidPred(:,1),'--g','Color',rand(1,3),'LineWidth',2);
    
    coll = [255,0,0;0,255,0;185,100,150;0,255,255;150,0,255;...
        255,0,255;255,150,0;200,90,150;90,140,90;0,0,255]./255;
    
    figure; hold on;
    count = 1;
    for j = 2.5:2.5:size(errT,3)
        errTH = errT(:,:,round(j));
        activeMeas = ~cellfun(@isempty,{errTH(:).errCentroidN});
        plot(xdim(activeMeas),[errTH(:).errCentroidN],'-x','Color',coll(count,:),'LineWidth',1.5,'MarkerSize',5);
        min_err(count) = min([errTH(:).errCentroidN]);
        max_err(count) = max([errTH(:).errCentroidN]);
        mean_err(count) = mean([errTH(:).errCentroidN]);
        med_err(count) = median([errTH(:).errCentroidN]);
        minute_c(count) = j*obj_r.frameRate*3/60;
        legendStr(count) = {strcat('Horizon :', num2str(j*obj_r.frameRate*3/60),' minutes.')};
        count = count + 1;
    end
    set(gca,'FontSize',11);
    set(gca,'YTick',[0:25:400]);
    title('Errors for multiple horizons (in pixels on the image plane).');
    xlabel('Frame Number for the Tracking Algorithm');
    ylabel('Error (in pixels on the image plane)');
    grid on;
    legend(legendStr,'Location','Northeast');
    
    figure, plot(minute_c,max_err,'-x',minute_c,mean_err,'-x',minute_c,min_err,'-x','LineWidth',2,'MarkerSize',5);
    set(gca,'FontSize',13);
    title('Statistical Error for different horizons');
    xlabel('Horizon (minutes)');
    ylabel('Error (pixels)');
    set(gca,'YTick',[0:25:400]);
    set(gca,'XTick',[0:0.5:5]);
    grid on;
    legend('Maximum','Mean','Minimum','Location','Northwest');
    
    
    figure; hold on;
    count = 1;
    for j = 2.5:2.5:size(errT,3)
        errTH = errT(:,:,round(j));
        activeMeas = ~cellfun(@isempty,{errTH(:).errCentroidProjN});
        plot(xdim(activeMeas),[errTH(:).errCentroidProjN],'-x','Color',coll(count,:),'LineWidth',1.5,'MarkerSize',5);
        min_errP(count) = min([errTH(:).errCentroidProjN]);
        max_errP(count) = max([errTH(:).errCentroidProjN]);
        mean_errP(count) = mean([errTH(:).errCentroidProjN]);
        med_errP(count) = median([errTH(:).errCentroidProjN]);
        minute_cP(count) = j*obj_r.frameRate*3/60;
        %legendStr(count) = {strcat('Prediction error for horizon of: ', num2str(j*obj_r.frameRate*3/60),' minutes.')};
        legendStr(count) = {strcat('Horizon :', num2str(j*obj_r.frameRate*3/60),' minutes.')};
        count = count + 1;
    end
    set(gca,'FontSize',11);
    set(gca,'YTick',[0:50:1400]);
    title('Errors for multiple horizons (in meters on the projected plane).');
    xlabel('Frame Number for the Tracking Algorithm');
    ylabel('Error (in meters on the projected plane)');
    grid on;
    legend(legendStr,'Location','Northwest');
    
    %,minute_cP,med_errP,'-x'
    figure, plot(minute_cP,max_errP,'-x',minute_cP,mean_errP,'-x',minute_cP,min_errP,'-x','LineWidth',2,'MarkerSize',5);
    set(gca,'FontSize',13);
    title('Statistical Error for different horizons');
    xlabel('Horizon (minutes)');
    ylabel('Error (meters)');
    set(gca,'YTick',[0:50:1400]);
    set(gca,'XTick',[0:0.5:5]);
    grid on;
    legend('Maximum','Mean','Minimum','Location','Northwest');
end


%%
ind_list = [94,231,224,266,272];
for ii = 1:length(ind_list)
    errT = errs(:,ind_list(ii),:);
    count = 1;
    for j = 2.5:2.5:size(errT,3)
        errTH = errT(:,:,round(j));
        %activeMeas = ~cellfun(@isempty,{errTH(:).errCentroidN});
        %plot(xdim(activeMeas),[errTH(:).errCentroidN],'-x','Color',coll(count,:),'LineWidth',1.5,'MarkerSize',5);
        min_err(ii,count) = min([errTH(:).errCentroidN]);
        max_err(ii,count) = max([errTH(:).errCentroidN]);
        mean_err(ii,count) = mean([errTH(:).errCentroidN]);
        med_err(ii,count) = median([errTH(:).errCentroidN]);
        minute_c(ii,count) = j*obj_r.frameRate*3/60;
        if(ind_list*(ii)==266)
            min_err(ii,count) = min([errTH(1:174).errCentroidN]);
            max_err(ii,count) = max([errTH(1:174).errCentroidN]);
            mean_err(ii,count) = mean([errTH(1:174).errCentroidN]);
            med_err(ii,count) = median([errTH(1:174).errCentroidN]);
        end
        %legendStr(count) = {strcat('Horizon :', num2str(j*obj_r.frameRate*3/60),' minutes.')};
        count = count + 1;
    end
    count = 1;
    for j = 2.5:2.5:size(errT,3)
        errTH = errT(:,:,round(j));
        %activeMeas = ~cellfun(@isempty,{errTH(:).errCentroidProjN});
        %plot(xdim(activeMeas),[errTH(:).errCentroidProjN],'-x','Color',coll(count,:),'LineWidth',1.5,'MarkerSize',5);
        min_errP(ii,count) = min([errTH(:).errCentroidProjN]);
        max_errP(ii,count) = max([errTH(:).errCentroidProjN]);
        mean_errP(ii,count) = mean([errTH(:).errCentroidProjN]);
        med_errP(ii,count) = median([errTH(:).errCentroidProjN]);
        minute_cP(ii,count) = j*obj_r.frameRate*3/60;
        if(ind_list*(ii)==266)
            min_err(ii,count) = min([errTH(1:174).errCentroidN]);
            max_err(ii,count) = max([errTH(1:174).errCentroidN]);
            mean_err(ii,count) = mean([errTH(1:174).errCentroidN]);
            med_err(ii,count) = median([errTH(1:174).errCentroidN]);
        end
        %legendStr(count) = {strcat('Prediction error for horizon of: ', num2str(j*obj_r.frameRate*3/60),' minutes.')};
        %legendStr(count) = {strcat('Horizon :', num2str(j*obj_r.frameRate*3/60),' minutes.')};
        count = count + 1;
    end
end

min_err_min = min(min_err);
mean_err_min = min(mean_err);
max_err_min = min(max_err);
min_err_max = max(min_err);
mean_err_max = max(mean_err);
max_err_max = max(max_err);
min_err_mean = mean(min_err);
mean_err_mean = mean(mean_err);
max_err_mean = mean(max_err);
minute_c = minute_c(1,:);

figure, plot(minute_c,max_err_min,'-x',minute_c,mean_err_min,'-x',minute_c,min_err_min,'-x','LineWidth',2,'MarkerSize',5);
set(gca,'FontSize',13);
title('Statistical Error for different horizons (Best)');
xlabel('Horizon (minutes)');
ylabel('Error (pixels)');
set(gca,'YTick',[0:5:400]);
set(gca,'XTick',[0:0.5:5]);
grid on;
legend('Maximum','Mean','Minimum','Location','Northwest');

figure, plot(minute_c,max_err_mean,'-x',minute_c,mean_err_mean,'-x',minute_c,min_err_mean,'-x','LineWidth',2,'MarkerSize',5);
set(gca,'FontSize',13);
title('Statistical Error for different horizons (Avg.)');
xlabel('Horizon (minutes)');
ylabel('Error (pixels)');
set(gca,'YTick',[0:15:400]);
set(gca,'XTick',[0:0.5:5]);
grid on;
legend('Maximum','Mean','Minimum','Location','Northwest');

min_errP_min = min(min_errP);
mean_errP_min = min(mean_errP);
max_errP_min = min(max_errP);
min_errP_max = max(min_errP);
mean_errP_max = max(mean_errP);
max_errP_max = max(max_errP);
min_errP_mean = mean(min_errP);
mean_errP_mean = mean(mean_errP);
max_errP_mean = mean(max_errP);
minute_cP = minute_cP(1,:);

figure, plot(minute_cP,max_errP_min,'-x',minute_cP,mean_errP_min,'-x',minute_cP,min_errP_min,'-x','LineWidth',2,'MarkerSize',5);
set(gca,'FontSize',13);
title('Statistical Error for different horizons (Best)');
xlabel('Horizon (minutes)');
ylabel('Error (meters)');
set(gca,'YTick',[0:30:1400]);
set(gca,'XTick',[0:0.5:5]);
grid on;
legend('Maximum','Mean','Minimum','Location','Northwest');

figure, plot(minute_cP,max_errP_mean,'-x',minute_cP,mean_errP_mean,'-x',minute_cP,min_errP_mean,'-x','LineWidth',2,'MarkerSize',5);
set(gca,'FontSize',13);
title('Statistical Error for different horizons (Avg.)');
xlabel('Horizon (minutes)');
ylabel('Error (meters)');
set(gca,'YTick',[0:50:1400]);
set(gca,'XTick',[0:0.5:5]);
grid on;
legend('Maximum','Mean','Minimum','Location','Northwest');