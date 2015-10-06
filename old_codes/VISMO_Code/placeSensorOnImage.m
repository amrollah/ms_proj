

load 'calib_model_01_2014_v2.mat';
load 'extrinsic_calib_model_01_2014.mat';
cbh = 1500;
%% generic parameters of model3d
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

frame = imread('2014_06_12_11_42_44_591.jpeg');
DN = datenum([2014,06,12,11,42,44]);
frame(~image_mask) = 0;  % we eliminate the parts that are angularly behind the camera observation point
res = cv.matchTemplate(frame,sunPattern); % find the position of the sun in the current frame through template matching with a previously obtained sun template
[y,x]=ind2sub(size(res),find(res==min(min(res)))); % find the pixel sun position on the image
s_s = size(sunPattern);
x_s = s_s(2);
y_s = s_s(1);
[X,Y] = meshgrid(1:ocam_model.height,1:ocam_model.width);
sunP = sqrt((X-(x+x_s/2)).^2+(Y-(y+y_s/2)).^2)<=180;
frame = cv.circle(frame,[x+x_s/2,y+y_s/2],80,'Color',[0,0,0],'Thickness',-1); % eliminate the sun putting a black circle on top of it.
r_b = double(frame(:,:,1))./double(frame(:,:,3)); % R/B thresholding applied
rgb = double(frame(:,:,1))+double(frame(:,:,2))+double(frame(:,:,3)); % to distinguish other objects like buildings from clouds, a test on brightness
mask = r_b > 0.72 & rgb > 100;             %% brightness auto thing is needed because of this ..
g_b = double(frame(:,:,2))./double(frame(:,:,3));
mask2 = g_b > 0.8;
mask(sunP) = mask2(sunP);
mask = medfilt2(mask,[10 10]);
contours = cv.findContours(mask); % contour extraction using openCV method
contoursProj = cellfun(@(x) num2cell(GetImage2Real(cell2mat(x(:)),ocam_model,ecam_model,cbh),2)',contours,'UniformOutput',false);
centroids = cellfun(@(x) (min(cell2mat(x(:)),[],1)+max(cell2mat(x(:)),[],1))/2,contours,'UniformOutput',false); % get the center of rectengular bounding box around the contour
centroidsProj = cellfun(@(x) (min(cell2mat(x(:)),[],1)+max(cell2mat(x(:)),[],1))/2,contoursProj,'UniformOutput',false);
bboxes = cellfun(@(x) [min(cell2mat(x(:)),[],1),max(cell2mat(x(:)),[],1)-min(cell2mat(x(:)),[],1)] ,contours,'UniformOutput',false); % form the bounding boxes using the extreme points
areas = cell2mat(cellfun(@(x) cv.contourArea(x(:)) ,contours,'UniformOutput',false))';
%areas = cell2mat(cellfun(@(x) prod(max(cell2mat(x(:)),[],1)-min(cell2mat(x(:)),[],1)) ,contours,'UniformOutput',false))';

% elminate the objects with contour smaller than 350.
consider = areas > 0; % if the area of an object is less than 350 pixels, we discard that

% take only the objects that we want to consider
centroids = centroids(consider);
bboxes = bboxes(consider);
areas = areas(consider);
contours = contours(consider);

centroids = vertcat(centroids{:});
centroidsProj = GetImage2Real(centroids(:,1:2),ocam_model,ecam_model,cbh); % the projected centroids are found using the mapping

bboxes = vertcat(bboxes{:}); % put all the boxes in a 2D array, so that the annotation algorithm can display them
bboxesTemp1 = GetImage2Real(bboxes(:,1:2),ocam_model,ecam_model,cbh); % upper corner of bboxes in projected world
bboxesTemp2 = GetImage2Real(bboxes(:,1:2)+bboxes(:,3:4),ocam_model,ecam_model,cbh); % lower corner of bboxes in projected world
bboxesTempDim = abs(bboxesTemp2 - bboxesTemp1); % dimensions of bboxes (in a sense)

Time = pvl_maketimestruct(DN, UTC);
[SunAz, ~, ApparentSunEl, ~] = pvl_ephemeris(Time, Location);
zenith = 90-ApparentSunEl;
azimuth = SunAz;
[x,y,~] = sph2cart((90-azimuth)*deg2radf,(90-zenith)*deg2radf,cbh/cos(zenith*deg2radf));
shadowDisp_d = [-x,-y];
        
frame2 = frame;
sensor_location1 = GetRelativePosition(Location,SensorLocation)-shadowDisp_d;
sensor_location1_on_image = GetReal2Image(sensor_location1,ocam_model,ecam_model,cbh);
frame2 = cv.circle(frame2,sensor_location1_on_image,10,'Color',[255,0,0],'Thickness',-1);
sensor_location2 = GetRelativePosition(Location,SensorLocation2)-shadowDisp_d;
sensor_location2_on_image = GetReal2Image(sensor_location2,ocam_model,ecam_model,cbh);
frame2 = cv.circle(frame2,sensor_location2_on_image,10,'Color',[0,0,255],'Thickness',-1);
sensor_location3 = GetRelativePosition(Location,SensorLocation3)-shadowDisp_d;
sensor_location3_on_image = GetReal2Image(sensor_location3,ocam_model,ecam_model,cbh);
frame2 = cv.circle(frame2,sensor_location3_on_image,10,'Color',[255,255,0],'Thickness',-1);
sensor_location4 = GetRelativePosition(Location,SensorLocation4)-shadowDisp_d;
sensor_location4_on_image = GetReal2Image(sensor_location4,ocam_model,ecam_model,cbh);
frame2 = cv.circle(frame2,sensor_location4_on_image,10,'Color',[255,0,255],'Thickness',-1);
sensor_location5 = GetRelativePosition(Location,SensorLocation5)-shadowDisp_d;
sensor_location5_on_image = GetReal2Image(sensor_location5,ocam_model,ecam_model,cbh);
frame2 = cv.circle(frame2,sensor_location5_on_image,10,'Color',[255,255,255],'Thickness',-1);
sensor_location6 = GetRelativePosition(Location,SensorLocation6)-shadowDisp_d;
sensor_location6_on_image = GetReal2Image(sensor_location6,ocam_model,ecam_model,cbh);
frame2 = cv.circle(frame2,sensor_location6_on_image,10,'Color',[0,0,0],'Thickness',-1);

figure, imshow(frame2);