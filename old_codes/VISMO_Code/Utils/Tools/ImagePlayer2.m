%% Video Player for specific use.
% The user can play-pause, fast forward-backward the video using the
% keyboard. Also, implementing a metric on the frame, this player can keep
% track of this metric. 
%
%Spesific purpose video player with variable frame-rate capability:
%
%Inputs:
% video_name: path to the input video file
% live: option indicating that a live stream is being played
% frameRate: option indication the number of frames that are going to be
% skipped for every displayed
%
%Outputs:
% ratio: a matrix with two columns, first is the date and the second is the
% result of the metric under consideration

% B. Zeydan, 19. Nov. 2013

function ImagePlayer2(image_folder,live,frameRate,ocam_model,ecam_model,data1,data2,data3,data4)
% create system objects used for 2reading video, detecting moving objects,
% and displaying the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Location.latitude =  47.459223;    %   47°27'33.92" N  Local latitude of CHCRC.C1 roof  (positive if north of equator) 47.459223 47+27/60+33.92/60/60
Location.longitude = 8.277023;     %   8°16'38.50" E   Local Longitude of CHCRC.C1 roof (negative if west) 8.277023 8+16/60+38.50/60/60
Location.altitude =  455;                     %   [m]
Location.distance = 0;                        %   [m] Distance to the camera
Location.bearing = 0;                         %   0° 0' 0" Bearing to the camera wrt the North
%sensor location 2 47.458385, 8.279722
SensorLocation.latitude =  47.459519;    %   47° 27' 34.27" N  Local latitude of CHCRC.C1 roof  (positive if north of equator) 47.459519 47+27/60+34.27/60/60
SensorLocation.longitude = 8.277495;     %   8° 16' 38.98" E    Local Longitude of CHCRC.C1 roof (negative if west) 8.277495 8+16/60+38.98/60/60
SensorLocation.altitude =  455;                     %   [m]
SensorLocation.distance = 48;                       %   [m] Distance to the camera
SensorLocation.bearing = 47+9/60+10/60/60;          %   47° 09' 10" Bearing to the camera wrt the North 47+9/60+10/60/60

SensorLocation2.latitude = 47.458385;
SensorLocation2.longitude = 8.279722;
SensorLocation2.altitude = 455;
SensorLocation2.distance = -1;
SensorLocation2.bearing = -1;

SensorLocation3.latitude = 47.458450;
SensorLocation3.longitude = 8.281255;
SensorLocation3.altitude = 455;
SensorLocation3.distance = -1;
SensorLocation3.bearing = -1;

%47.481 8.293
SensorLocation4.latitude = 47.481;
SensorLocation4.longitude = 8.293;
SensorLocation4.altitude = 455;
SensorLocation4.distance = -1;
SensorLocation4.bearing = -1;

%47.411 8.316
SensorLocation5.latitude = 47.411;
SensorLocation5.longitude = 8.316;
SensorLocation5.altitude = 455;
SensorLocation5.distance = -1;
SensorLocation5.bearing = -1;

%47.505 8.248
SensorLocation6.latitude = 47.505;
SensorLocation6.longitude = 8.248;
SensorLocation6.altitude = 455;
SensorLocation6.distance = -1;
SensorLocation6.bearing = -1;

UTC = 2;
deg2radf = pi/180;
cbh = 1500;
% here we set-up the relevant objects.
obj = setupSystemObjects(image_folder,live,frameRate,data1,data2,data3,data4);

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
        frame2 = frame;
        
        Time = pvl_maketimestruct(obj.frameDatenum, UTC);
        
        [SunAz, ~, ApparentSunEl, ~] = pvl_ephemeris(Time, SensorLocation);
        zenith = 90-ApparentSunEl;
        azimuth = SunAz;
        [x,y,~] = sph2cart((90-azimuth)*deg2radf,(90-zenith)*deg2radf,cbh/cos(zenith*deg2radf));
        shadowDisp_d1 = [-x,-y];
        sensor_location1 = GetRelativePosition(Location,SensorLocation)-shadowDisp_d1;
        sensor_location1_on_image = GetReal2Image(sensor_location1,ocam_model,ecam_model,cbh);
        frame2 = cv.circle(frame2,sensor_location1_on_image,15,'Color',[255,0,0],'Thickness',-1);
        
        [SunAz, ~, ApparentSunEl, ~] = pvl_ephemeris(Time, SensorLocation2);
        zenith = 90-ApparentSunEl;
        azimuth = SunAz;
        [x,y,~] = sph2cart((90-azimuth)*deg2radf,(90-zenith)*deg2radf,cbh/cos(zenith*deg2radf));
        shadowDisp_d2 = [-x,-y];
        sensor_location2 = GetRelativePosition(Location,SensorLocation2)-shadowDisp_d2;
        sensor_location2_on_image = GetReal2Image(sensor_location2,ocam_model,ecam_model,cbh);
        frame2 = cv.circle(frame2,sensor_location2_on_image,15,'Color',[0,0,255],'Thickness',-1);
        frame2 = cv.putText(frame2,'S3',sensor_location2_on_image-[20,20],'Color',[0,0,255],'Thickness',3);
        
        [SunAz, ~, ApparentSunEl, ~] = pvl_ephemeris(Time, SensorLocation3);
        zenith = 90-ApparentSunEl;
        azimuth = SunAz;
        [x,y,~] = sph2cart((90-azimuth)*deg2radf,(90-zenith)*deg2radf,cbh/cos(zenith*deg2radf));
        shadowDisp_d3 = [-x,-y];
        sensor_location3 = GetRelativePosition(Location,SensorLocation3)-shadowDisp_d3;
        sensor_location3_on_image = GetReal2Image(sensor_location3,ocam_model,ecam_model,cbh);
        frame2 = cv.circle(frame2,sensor_location3_on_image,15,'Color',[255,255,0],'Thickness',-1);
        
        [SunAz, ~, ApparentSunEl, ~] = pvl_ephemeris(Time, SensorLocation4);
        zenith = 90-ApparentSunEl;
        azimuth = SunAz;
        [x,y,~] = sph2cart((90-azimuth)*deg2radf,(90-zenith)*deg2radf,cbh/cos(zenith*deg2radf));
        shadowDisp_d4 = [-x,-y];
        sensor_location4 = GetRelativePosition(Location,SensorLocation4)-shadowDisp_d4;
        sensor_location4_on_image = GetReal2Image(sensor_location4,ocam_model,ecam_model,cbh);
        frame2 = cv.circle(frame2,sensor_location4_on_image,15,'Color',[255,0,255],'Thickness',-1);
        
        [SunAz, ~, ApparentSunEl, ~] = pvl_ephemeris(Time, SensorLocation5);
        zenith = 90-ApparentSunEl;
        azimuth = SunAz;
        [x,y,~] = sph2cart((90-azimuth)*deg2radf,(90-zenith)*deg2radf,cbh/cos(zenith*deg2radf));
        shadowDisp_d5 = [-x,-y];
        sensor_location5 = GetRelativePosition(Location,SensorLocation5)-shadowDisp_d5;
        sensor_location5_on_image = GetReal2Image(sensor_location5,ocam_model,ecam_model,cbh);
        frame2 = cv.circle(frame2,sensor_location5_on_image,15,'Color',[255,255,255],'Thickness',-1);
        frame2 = cv.putText(frame2,'S2',sensor_location5_on_image-[20,20],'Color',[255,255,255],'Thickness',3);
        
        [SunAz, ~, ApparentSunEl, ~] = pvl_ephemeris(Time, SensorLocation6);
        zenith = 90-ApparentSunEl;
        azimuth = SunAz;
        [x,y,~] = sph2cart((90-azimuth)*deg2radf,(90-zenith)*deg2radf,cbh/cos(zenith*deg2radf));
        shadowDisp_d6 = [-x,-y];
        sensor_location6 = GetRelativePosition(Location,SensorLocation6)-shadowDisp_d6;
        sensor_location6_on_image = GetReal2Image(sensor_location6,ocam_model,ecam_model,cbh);
        frame2 = cv.circle(frame2,sensor_location6_on_image,15,'Color',[255,255,255],'Thickness',-1);
        frame2 = cv.putText(frame2,'S1',sensor_location6_on_image-[20,20],'Color',[255,255,255],'Thickness',3);
        
        frame = frame2;
        displayTrackingResults(); % here we display the status and results
        obj.request = 0;
    else
        displayTrackingResults();
    end
    ratio = obj.ratio;
end

    function obj = setupSystemObjects(image_folder,live,frameRate,data1,data2,data3,data4)
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
        obj.plotPlayer = ModPlayer([740, 500, 700, 300],'plot'); % for the GHI-value plot
        obj.frameRate = frameRate;
        obj.paused = 0;
        obj.request = 0;
        obj.ratio = [];
        obj.frameDateNum = 0;
        obj.data1 = data1;
        obj.data2 = data2;
        obj.data3 = data3;
        obj.data4 = data4;
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
                frame = imread(char(strcat(obj.reader, '/', obj.im_list(obj.frameCount))));  % read the image directly, no need for the lock operations
                obj.frameCount = obj.frameCount + obj.frameRate;
            else
                frame = []; % if no more frames in the recording, directly terminate
                display('Done reading strored images... Terminating...');
            end
        end
        
        if(~isempty(frame))
            r_b = double(frame(:,:,3))-double(frame(:,:,1)); % R/B thresholding applied
            obj.ratio(end+1,:) = [obj.frameDatenum, sum(r_b(:))/1000];
        end
        minuteLag = 15; % go back for this much of minutes
        queryDates = datestr([-minuteLag:1:minuteLag]'/(60*24)+obj.frameDatenum);
        irradianceDB = resample(obj.data1,queryDates);
        obj.irradianceValues1 = irradianceDB.Data(1:end);
        obj.irradianceTimes1 = [-minuteLag:1:minuteLag]';%obj.occlusionTimes';
        irradianceDB2 = resample(obj.data2,queryDates);
        obj.irradianceValues2 = irradianceDB2.Data(1:end);
        obj.irradianceTimes2 = [-minuteLag:1:minuteLag]';%obj.occlusionTimes';
        irradianceDB3 = resample(obj.data3,queryDates);
        obj.irradianceValues3 = irradianceDB3.Data(1:end);
        obj.irradianceTimes3 = [-minuteLag:1:minuteLag]';%obj.occlusionTimes';
        if(isempty(obj.irradianceValues1)) % if not found the associated entries, set to zero
            obj.irradianceValues1 = 0;
            obj.irradianceTimes1 = 0;
            obj.irradianceValues2 = 0;
            obj.irradianceTimes2 = 0;
            obj.irradianceValues3 = 0;
            obj.irradianceTimes3 = 0;
        end
        cbhDB = resample(obj.data4, datestr(obj.frameDatenum)); % query the values
        cbh =  cbhDB.Data(1);
    end

    function displayTrackingResults()
        obj.videoPlayer.Step(frame,char(strcat(obj.frameName,' Frame Rate: ',num2str(obj.frameRate))));
        %obj.plotPlayer.Step(obj.ratio,strcat('Time Relative to ', char(obj.frameName), ' (Min.)') ,'Irradiance (W/m2)','Irradiance Value Changes');
        obj.plotPlayer.Step([obj.irradianceTimes1,obj.irradianceValues1,obj.irradianceTimes2,obj.irradianceValues2,obj.irradianceTimes3,obj.irradianceValues3],strcat('Time Relative to ', datestr(obj.frameDatenum), ' (Min.)') ,'Irradiance (W/m2)','Irradiance Values',['S1 ';'S2 ';'S3 ']);
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
                        obj.frameCount = obj.frameCount+obj.frameRate;
                        obj.request = 1;
                    end
                case 'downarrow'
                    if(obj.frameRate>1)
                        obj.frameRate = obj.frameRate-1;
                    end
                case 'l'
                    obj.paused = ~obj.paused;
                otherwise
            end
        end
    end
end