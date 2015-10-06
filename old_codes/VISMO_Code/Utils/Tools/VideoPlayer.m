%% Video Player for generic use.
% The user can play-pause, fast forward-backward the video using the
% keyboard.
% Generic purpose video player with variable frame-rate capability:
%Inputs:
% video_name: path to the input video file
% live: option indicating that a live stream is being played
% frameRate: option indication the number of frames that are going to be
% skipped for every displayed

%Outputs:
% NONE

% B. Zeydan, 28. Sep. 2013

function VideoPlayer(image_folder,live,frameRate)
% create system objects used for 2reading video, detecting moving objects,
% and displaying the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cbh = 300;
% here we set-up the relevant objects.
obj = setupSystemObjects(image_folder,live,frameRate);

frame = [];
if(isempty(obj.reader))
    disp('Camera does not respond. Terminating...');
    return;
end

% This is the main loop of the program
while obj.frameCount < obj.finalFrame
    if(~obj.paused  || obj.request)
        frame = readFrame();  % readFrame function returns the image in order from the buffer folder or the recording folder that is determined during the setup
        if(isempty(frame))    % if the folder does not exist or the frame is not correct, we immediately terminate
            disp('Camera does not respond. Terminating...');
            close all;
            return;
        end
        displayTrackingResults(); % here we display the status and results
        obj.request = 0;
    else
        displayTrackingResults();
    end
end

    function obj = setupSystemObjects(image_folder,live,frameRate)
        % Initialize Video I/O
        % Create objects for reading a video from a file, drawing the tracked
        % objects in each frame, and playing the video.
        
        % create a video file reader
        try_ct = 0;
        obj.live = live;
        
        % if the read functionality is live, than we first check if the
        % buffer folder exists so the camera is already operational
        % if it cannot find the buffer folder, then it runs a bussy-wait
        % procedure on it. When it times out, it terminates the setup
        % function and terminates the program before ever proceeding. if
        % the read functionality is not live, then is is directly
        % initiated to be checked during the first read.
        if(obj.live)
            while(1)
                try_ct = try_ct + 1;
                if(exist(char([image_folder, '/Buffer']))==7)
                    obj.reader = char([image_folder, '/Buffer']);
                    break;
                else
                    if(try_ct > 10)
                        obj.reader = [];
                        disp('No input stream found.. Exit.');
                        return;
                    else
                        disp(char(['Trying to read for: ' num2str(try_ct)]));
                        pause(2);
                    end
                end
            end
            obj.vid_ims = [];
            obj.im_list = [];
            obj.finalFrame = inf;
        else
            obj.reader = image_folder;
            obj.vid_ims = ls(obj.reader); % get the list of files in the folder
            obj.im_list = cellstr(obj.vid_ims(3:end,:)); % get the name list of those list of files
            obj.finalFrame = length(obj.im_list);
        end
        
        
        obj.frameCount = 1; % the order of frames that are processed.
        obj.frameName = []; % the file name of the frame that is processed the last.
        obj.videoPlayer = VidPlayer([20, 500, 700, 400],@cb); % for the real image
        obj.frameRate = frameRate;
        obj.paused = 0;
        obj.request = 0;
        obj.recording = 0;
        obj.vid_id = 0;
        obj.vid_n = [];
        obj.vidObj = [];
        obj.frameID = [];
        obj.frameDatenum = 0;
    end

    function frame = readFrame()
        % in this function we read the frame on the order from the folder
        
        % if the destination buffer do not exist, than we immediately
        % terminate
        if(exist(obj.reader) ~=7)
            disp('Input stream should have been closed... Terminating...');
            frame = []; % indicator of an unsuccessful read trial
            return;
        end
        
        
        read_try = 1; % try to read a file for the first time
        if(obj.live) % if this is the live mode, meaning, we are going at the same paste as the camera
            vid_ims = ls(obj.reader); % get the list of files in the folder
            im_list = cellstr(vid_ims(3:end,:)); % get the name list of those list of files
            while(length(im_list{1}) == 0) % wait for an image to be written to the buffer
                if(exist(obj.reader) ~=7) % if the buffer folder is deleted during the wait, this means the camera has been shut down, so terminate
                    disp('Input stream should have been closed... Terminating...');
                    frame = []; % indicator of an unsuccessful read trial
                    return;
                end
                vid_ims = ls(obj.reader); % get the list of files in the folder
                im_list = cellstr(vid_ims(3:end,:)); % get the name list of those list of files
                read_try = read_try+1; % increment the number of trials and go for another try
                disp(char(['Trying to read for: ' num2str(read_try)])); % report the action to the user.
                if(read_try > 8) % if a read is tried for more than 8 times, then terminate with an unsuccessfull read
                    disp('Input stream not available... Terminating...');
                    frame = []; % unsuccessful read indicator
                    return;
                end
                pause(1); % wait 1 msec between two read trials, we can increase this time unit, but not so necessary
            end
            % wait until the name of the file in the buffer folder is
            % different from the previous read to make sure that the
            % observation frame is updated
            while(strcmp(char(strcat(obj.reader, '/', im_list(1))),obj.frameName)==1)
                pause(1);
                vid_ims = ls(obj.reader);
                im_list = cellstr(vid_ims(3:end,:));
            end
            % if the code is here, this means that the new file is ready
            % to be read
            fid = fopen(char(strcat(obj.reader, '/', im_list(1))), 'r'); % we try to read in case there is an operating system related lock on the file. This may happen if the camera is writing the
            % the file as we try to read it. if this application manages to open the file, than the file is locked by it and can be read.
            while (fid == -1) % if the camera is writing the file then apply the same trial based procedure until the write is over.
                if(read_try > 8)
                    disp('Input stream not available... Terminating...');
                    frame = [];
                    return;
                end
                disp(['Retrying from this one...   ' char(strcat(obj.reader, '/', im_list(1)))]);
                vid_ims = ls(obj.reader);
                im_list = cellstr(vid_ims(3:end,:));
                fid = fopen(char(strcat(obj.reader, '/', im_list(1))), 'r'); % try to lock the file once again
                read_try = read_try + 1;
                pause(0.5);
            end
            
            frame = imread(char(strcat(obj.reader, '/', im_list(1)))); % since the file is locked, it is safe to read it.
            fclose(fid); % since the file is read, the lock can be released
            
            obj.frameName = char(strcat(obj.reader, '/', im_list(1))); % record the name of the file for the check operation on the up-to-date information of the file
            obj.frameCount = obj.frameCount + 1; % increase the frame count
        else
            % if it is not a live operation, directly read the file in
            % order
            if(length(obj.im_list) >= obj.frameCount) % if there are still some images to be read
                date_vec = strread(char(obj.im_list(obj.frameCount)),'%s','delimiter','_');
                obj.frameDatenum = datenum(str2num(char(date_vec(1))),str2num(char(date_vec(2))),str2num(char(date_vec(3))),str2num(char(date_vec(4))),str2num(char(date_vec(5))),str2num(char(date_vec(6))));
                obj.frameName = strcat('Date:  ',date_vec(3),'.',date_vec(2),'.',date_vec(1),' Hour: ',date_vec(4),':',date_vec(5),':',date_vec(6));
                obj.frameID = strcat(date_vec(3),'.',date_vec(2),'.',date_vec(1),' ,',date_vec(4),':',date_vec(5),':',date_vec(6));
                frame = imread(char(strcat(obj.reader, '/', obj.im_list(obj.frameCount))));  % read the image directly, no need for the lock operations
                if(obj.recording)
                    frame = cv.putText(frame,'R',[1700,100],'FontScale',2,'Thickness',7,'Color',[0,255,0]);
                end
                %frame = cv.putText(frame,char(obj.frameID),[10,1790],'FontScale',1.3,'Thickness',4,'Color',[0,255,0]);
                frame = cv.putText(frame,char(datestr(obj.frameDatenum)),[5,1790],'FontScale',1.2,'Thickness',3,'Color',[255,255,255]);
                obj.frameCount = obj.frameCount + obj.frameRate;
            else
                if(obj.recording)
                    close(obj.vidObj);
                end
                frame = []; % if no more frames in the recording, directly terminate
                display('Done reading strored images... Terminating...');
            end
        end
        if(obj.recording)
            writeVideo(obj.vidObj, frame);
        end
    end
    function displayTrackingResults()
        obj.videoPlayer.Step(frame,char(strcat(obj.frameName,' Frame Rate: ',num2str(obj.frameRate))));
    end
    function cb(~,evnt)
        if strcmp(evnt.Modifier{1},'control')
            switch evnt.Key
                case 'uparrow'
                    if(obj.frameCount+obj.frameRate+1<=obj.finalFrame)
                        obj.frameRate = obj.frameRate+1;
                    end
                case 'leftarrow'
                    if(obj.frameCount-2*obj.frameRate>=1)
                        obj.frameCount = obj.frameCount-2*obj.frameRate;
                        obj.request = 1;
                    end
                case 'rightarrow'
                    if(obj.frameCount+obj.frameRate<=obj.finalFrame)
                        obj.request = 1;
                    end
                case 'downarrow'
                    if(obj.frameRate>1)
                        obj.frameRate = obj.frameRate-1;
                    end
                case 'l'
                    obj.paused = ~obj.paused;
                case 'r'
                    if(~obj.recording)
                        obj.vid_id = obj.vid_id+1;
                        obj.vid_n = char(strcat('vid_', num2str(obj.vid_id), '.avi'));
                        obj.vidObj = VideoWriter(obj.vid_n);
                        obj.vidObj.Quality = 100;
                        obj.vidObj.FrameRate = 4;
                        open(obj.vidObj);
                    else
                        close(obj.vidObj);
                    end
                    obj.recording = ~obj.recording;
                case 'i'
                    obj.recording = ~obj.recording;
                otherwise
            end
        end
    end
end