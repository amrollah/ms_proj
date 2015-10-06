%% Script for collecting tracking results over the whole procedure
% this script calls the method with camera inputs and the calibration model
% resulting tracks are saved on the hard drive depending on the state of 
% the save_files variable. If called with load_files set, this loads the last 
% saved state and if kalman_test is set at the same time it runs the tracking 
% algorithm without reading the images. This functionality is to gain some time
% in repeated testing.
% B. Zeydan, 31. Mar. 2013

%%
save_files = 1;
load_files = 0;
kalman_test = 0;
video_name = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/05/14'; % the name of the stream file (video or a collection of images)
if(~exist('ocam_model'))
    display('No calibration model. Terminating...');
    return;
end
live = 0;
%%
if(load_files)
    e_val = exist('recordings');
    if ~(e_val == 7)
        error('no recording exists');
    end
    cd recordings;
    c_list = sort(cell2mat(cellfun(@(x) str2num(x),cellstr(ls),'UniformOutput',false)));
    
    if ~isempty(c_list)
        b_num = c_list(end);
    else
        error('no recording exists');
    end
    
    b_str = num2str(b_num);
    cd(b_str);
    e_val1 = exist('recording_obs');
    e_val2 = exist('recording_pred');
    e_val3 = exist('db_objs');
    e_val4 = exist('db_ass');
    if((e_val1 == 0) || (e_val2 == 0) || (e_val3 == 0) || (e_val4 == 0))
        e_val1 = exist('observed.mat');
        e_val2 = exist('predicted.mat');
        e_val3 = exist('observations.mat');
        e_val4 = exist('assignments.mat');
        if((e_val1 ~= 2) || (e_val2 ~= 2) || (e_val3 ~= 2) || (e_val4 ~= 2))
            error('incomplete recording please delete the latest unsuccessful recording');
        else
            load('observed.mat');
            load('predicted.mat');
            load('observations.mat');
            load('assignments.mat');
        end
    else
        warning('Already existing recordings. Not loading. To reload delete recording_obs, recording_preds variables from your workspace and try again.');
    end
    cd ../..
    if(kalman_test)
        [recording_obs,recording_pred] = multiObjectTrackingCV_Cam(video_name,1,1,5000,ocam_model,live);
    end
else
    [recording_obs,recording_pred,db_objs,db_ass] = multiObjectTrackingCV_Cam(video_name,1,0,2000,ocam_model,live);
    %% saving to file part...
    if(save_files)
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
        save('observed', 'recording_obs');
        save('predicted', 'recording_pred');
        save('observations', 'db_objs');
        save('assignments', 'db_ass');
        
        cd ../..
        
    end
end