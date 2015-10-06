% --------------------------------------------
%  3D MAP reconstruction for VISMO
%
% Note: You need the PVLIB library to run the code!
%
% 18.09.2013
% Luis Dominguez
% Control & Optimization
% ABB Corporate Research
% --------------------------------------------

%close all
%clf
%clc
%clear all
%load picture

% %Conversion to real coordinates using scaled measurements
% L1 = 28.7/1.9*50;  % m
% L2 = 15.3/1.9*50;  % m
% %L1 = 40.9/2*500;  % m
% %L2 = 21.1/2*500;  % m
% 
% 
% % % PV Center Coordinates on map
% %PV_xc = 15.5/1.9*50  %m
% %PV_yc = 8.1/1.9*50  %m
% 
% % % Remaping cordinates, PV as origin
% % L1_lb =  -PV_xc
% % L1_ub = L1-PV_xc
% %
% % L2_lb = -PV_yc
% % L2_ub = L2-PV_yc
% % %-L2_lb + L2_ub %Check
% 
% % Camera Coordinates on map
% %Cam_xc = 20.5/2*500;  %m
% %Cam_yc = 10.1/2*500;  %m
% Cam_xc = 14.5/1.9*50;  %m
% Cam_yc = 7.2 /1.9*50;  %m
% 
% % Remaping cordinates,  camera as origin
% x_lb =    - Cam_xc;
% x_ub = L1 - Cam_xc;
% 
% y_lb =    - Cam_yc;
% y_ub = L2 - Cam_yc;
% 
% % Read Image
% %CHCRC = imread('CHCRC_Bigger.jpg');
% CHCRC = imread('CHCRC.bmp');
L1 = 40.9/1.8*1000;  % m
L2 = 21.1/1.8*1000;  % m


% % PV Center Coordinates on map
% PV_xc = 15.5/1.9*50  %m
% PV_yc = 8.1/1.9*50  %m

% % Remaping cordinates, PV as origin
% L1_lb =  -PV_xc
% L1_ub = L1-PV_xc
%
% L2_lb = -PV_yc
% L2_ub = L2-PV_yc
% %-L2_lb + L2_ub %Check

% Camera Coordinates on map
Cam_xc = 20.5/1.8*1000;  %m
Cam_yc = 10.1/1.8*1000;  %m

% Remaping cordinates,  camera as origin
x_lb =    - Cam_xc;
x_ub = L1 - Cam_xc;

y_lb =    - Cam_yc;
y_ub = L2 - Cam_yc;

% Read Image
CHCRC = imread('CHCRC_Bigger2.jpg');
%figure(1)
%imagesc(CHCRC,'XData',[x_lb x_ub ],'YData',[y_lb y_ub  ]);
%imagesc(CHCRC,'XData',[x_lb x_ub ],'YData',[-y_ub  -y_lb ]);

figure(2)
s = surface('XData',[x_lb x_ub; x_lb x_ub],...
    'YData',[y_lb y_lb; y_ub y_ub],...
    'ZData',[0 0; 0 0],'CData',flipdim(CHCRC,1),...
    'FaceColor','texturemap','EdgeColor','none');

%text(5,-200,10,'South','HorizontalAlignment','left', 'BackgroundColor','w')
%text(5,180,10,'North','HorizontalAlignment','left', 'BackgroundColor','w')
%text(-375,10,10,'West','HorizontalAlignment','left', 'BackgroundColor','w')
%text(355,10,10,'East','HorizontalAlignment','left', 'BackgroundColor','w')

xlabel ( '--X axis--' );
ylabel ( '--Y axis--' );
zlabel ( '--Z axis--' );

% Plot coordinate line on 3D map
hold on, plot([0 0],[y_lb y_ub],'m')
hold on, plot([x_lb x_ub],[0 0],'m')
%hold on, scatter3(model3D.Sen_xs,model3D.Sen_ys,0)

[ x, y, z ] = sphere ( 50 );

scale_factor  = 200; %Cam_yc;

x2 = x*scale_factor;
y2 = y*scale_factor;
z2 = z*scale_factor;

z2(z2<0)=0; %form hemisphere

c = 0.20 * ones ( size ( z2 ) );

x3 = x*scale_factor*3;
y3 = y*scale_factor*3;
z3 = z*scale_factor*3;
z3(z3<0)=0; %form hemisphere
c3 = 0.40 * ones ( size ( z3 ) );

x4 = x./z*(scale_factor*3);
y4 = y./z*(scale_factor*3);
z4 = z./z*(scale_factor*3);
z4(z4<0)=0; %form hemisphere
c4 = 0.80 * ones ( size ( z4 ) );
%Plot hemispheric view of fish eye
esph = surf ( x2, y2, z2, c, 'EdgeColor', 'none');
alpha(esph,0.1)


%esph3 = surf ( x3, y3, z3, c3, 'EdgeColor', 'none');
%alpha(esph3,0.1)

%esph4 = surf ( x4, y4, z4, c4, 'EdgeColor', 'none');
%alpha(esph4,0.1)

view(-20,49)

% Generate data for sun location
Location.latitude =  47+27/60+33.15/60/60;    %   47°27'33.92" N  Local latitude of CHCRC.C1 roof  (positive if north of equator) for the camera +47° 27' 33.16"
Location.longitude = 8+16/60+37.30/60/60;     %   8°16'38.50" E   Local Longitude of CHCRC.C1 roof (negative if west) for the camera +8° 16' 37.28"
                                              %   also +47° 27' 33.15", +8° 16' 37.30" is the camera location
Location.altitude =  455;                     %   [m]  

%start_date = datenum(2013,07,15);
%start_date = datenum(2013,07,13);
start_date = datenum(2013,07,12);

%dt = 18/24+ 10/60/24;
%dt = 14/24+ 32/60/24; %2:32
%dt = 9/24
dt = 12/24;
%dt = 19/24;
%dt= 9/24: 1/24 : 19/24;

% Data is given in 10-min time series
DN = dt +  start_date;
nDN= numel(DN);

disp('Displaying sun trajectory for the following times:' );
disp(' ');
disp(datestr(DN,'mmmm dd, yyyy HH:MM'));

% Siwtzerland has a Central European Summer Time (CEST) +0200 UTC
UTC = 2; %
Time = pvl_maketimestruct(DN, UTC);

% Pre-calculations & conversion factors
rad2degf = 180/pi;
deg2radf = pi/180;

%Solar position calculations (azimuth, zenith and elevation)
[SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(Time, Location);
ApparentZenith = 90-ApparentSunEl;

%         figure(),  plot(dt,SunEl), datetick('x')
%         plot(dt,ApparentSunEl)
%         figure(), plot(dt,SolarTime), datetick('x')
%         plot(dt,ApparentZenith)
%         figure(3), plot(dt,SunAz), datetick('x')


zenith  = ApparentZenith';

%  -------- Notes on the importance of checking the coordinate system --------------------
% There are several conventions for the solar azimuth, however it is traditionally defined as the angle between a line due south and the shadow cast by a vertical rod on Earth. This convention states the angle is positive if the line is east of south and negative if it is west of south.[1][2] For example due east would be 90° and due west would be -90°. Another convention is the reverse; it also has the origin at due south, but measures angles clockwise, so that due east is now negative and west now positive.[3]
% However, despite tradition, the most commonly accepted convention for analyzing solar radiation, e.g. for solar energy applications, is clockwise from due north, so east is 90°, south is 180° and west is 270°. This is the definition used by NREL in their solar position calculators[4] and is also the convention used in the formulas presented here.

% In land navigation, azimuth is usually denoted alpha, and defined as a horizontal angle measured clockwise from a north base line or meridian.[3][4] Azimuth has also been more generally defined as a horizontal angle measured clockwise from any fixed reference plane or easily established base direction line.[5][6][7]
% Today the reference plane for an azimuth is typically true north, measured as a 0° azimuth, though other angular units (grad, mil) can be used. Moving clockwise on a 360 degree circle, east has azimuth 90°, south 180°, and west 270°. There are exceptions: some navigation systems use south as the reference plane. Any direction can be the plane of reference, as long as it is clearly defined.
% Quite commonly, azimuths or compass bearings are stated in a system in which either north or south can be the zero, and the angle may be measured clockwise or anticlockwise from the zero. For example, a bearing might be described as "(from) south, (turn) thirty degrees (toward the) east" (the words in brackets are usually omitted), abbreviated "S30°E", which is the bearing 30 degrees in the eastward direction from south, i.e. the bearing 150 degrees clockwise from north. The reference direction, stated first, is always north or south, and the turning direction, stated last, is east or west. The directions are chosen so that the angle, stated between them, is positive, between zero and 90 degrees. If the bearing happens to be exactly in the direction of one of the cardinal points, a different notation, e.g. "due east", is used instead.

%A number of different spherical coordinate systems

% *following other conventions are used outside mathematics. In a geographical coordinate system positions are measured in latitude, longitude and height or altitude. There are a number of different celestial coordinate systems based on different fundamental planes and with different terms for the various coordinates. The spherical coordinate systems used in mathematics normally use radians rather than degrees and measure the azimuthal angle counter-clockwise rather than clockwise.*
%----------------------------------------------------------------------------------------------

azimuth = (360 - SunAz') + 90; % This is a fake azimuth used to display the azitmuthal projection (and the computation of the sun coordininates in x-y-z) using the espherical coordinates equations. Note that in this eequations the zero-angle is EAST and the angle increases counterclockwise. On the other hand, SunAz is the azimuth considering NORTH as the zero-angle increasing clock-wise


r = 300; % Assume a given (radio) distance (towards sun) to calculate coordinate points

% Using sphereical coordinates
x_sun = r * sin(zenith*deg2radf) .* cos(azimuth*deg2radf);
y_sun = r * sin(zenith*deg2radf) .* sin(azimuth*deg2radf);
z_sun = r * cos(zenith*deg2radf);

% if azimuth < 180 & azimuth > 90
%     x_sun = r * sin(zenith*deg2radf) .* sin( (180 - azimuth) * deg2radf);
%     y_sun = r * sin(zenith*deg2radf) .* cos( azimuth * deg2radf);
%     z_sun = r * cos(zenith*deg2radf);
% end
%
% if azimuth < 270 & azimuth > 180
% x_sun = r * sin(zenith*deg2radf) .* sin( (azimuth - 180) * deg2radf);
% y_sun = r * sin(zenith*deg2radf) .* cos( (azimuth - 180) * deg2radf);
% z_sun = r * cos(zenith*deg2radf);
% end

DomeAltitude = 200; %m

%  Draw the focus, center, and pole points in black.
hold on,  scatter3 ( [0, 0], [0, 0], [0, DomeAltitude], 'filled', 'k', 'SizeData', 100 );

%  Draw a line through the focus (also center), and noth pole.
hold on,  line ( [0,0], [0,0], [0, DomeAltitude], 'Color', 'k', 'LineWidth', 2 )

%  Draw the image points in red.
hold on, scatter3 ( x_sun', y_sun', z_sun', 'filled', 'r' );

%  Draw a line indicating the incidence vector towards the focus.
hold on,  line ( [zeros(1,nDN);  x_sun], ...
    [zeros(1,nDN);  y_sun], ...
    [zeros(1,nDN);  z_sun], 'Color', 'r', 'LineWidth', 1 )

%  Draw a line of the projection of the sun on ground.
%here the z-component is zero
hold on,  line ( [zeros(1,nDN);  x_sun], ...
    [zeros(1,nDN);  y_sun], ...
    zeros(2,nDN), 'LineStyle','--','Color', 'r', 'LineWidth', 1)

%  Draw a line of the sun trajectory -here the z-component is zero
hold on,  line ( x_sun, ...
    y_sun, ...
    z_sun, 'LineStyle','-','Color', 'y', 'LineWidth', 2 );

title ( 'Stereographic Projection of Sun Trajectory' )



% ----------------------- Shadow Detection ----------------------------
L1 = 40.9/1.8*1000;  % m
L2 = 21.1/1.8*1000;  % m


% % PV Center Coordinates on map
% PV_xc = 15.5/1.9*50  %m
% PV_yc = 8.1/1.9*50  %m

% % Remaping cordinates, PV as origin
% L1_lb =  -PV_xc
% L1_ub = L1-PV_xc
%
% L2_lb = -PV_yc
% L2_ub = L2-PV_yc
% %-L2_lb + L2_ub %Check

% Camera Coordinates on map
Cam_xc = 20.5/1.8*1000;  %m
Cam_yc = 10.1/1.8*1000;  %m

% Remaping cordinates,  camera as origin
x_lb =    - Cam_xc;
x_ub = L1 - Cam_xc;

y_lb =    - Cam_yc;
y_ub = L2 - Cam_yc;

% Read Image
CHCRC = imread('CHCRC_Bigger2.jpg');
%load 'calib_modelGood.mat';
%load 'extrinsic_calib_model4.mat';
%load 'image_mask.mat';
load 'calib_model_01_2014_v2.mat';
load 'extrinsic_calib_model_01_2014.mat';

image_mask = ones(ocam_model.height,ocam_model.width);
[X,Y] = meshgrid(1:ocam_model.height,1:ocam_model.width);
obj.X = X;
obj.Y = Y;
P = [X(:)';Y(:)'];
p = cam2world(P,ocam_model);
p = ecam_model.R*p;
maskTemp = p(3,:)<-0.15;
indicesTemp = sub2ind([ocam_model.height,ocam_model.width],P(1,:),P(2,:));
image_mask(indicesTemp) = maskTemp;
image_mask = logical(repmat(image_mask,[1,1,3]));
sunPattern = imread('sunPattern.jpg');

% Assume we know the position of the Cloud
cbh = 500;  % provided by the calibration model, cbh is the assumed height
% Here the cbh has been scaled down for visual inspection purpuses
frame = imread('oimage5.jpeg');
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

zenith_cloud = zenith(1);
azimuth_cloud = SunAz(1)';


Time = pvl_maketimestruct(DN, model3D.UTC);
[SunAz, ~, ApparentSunEl, ~]=pvl_ephemeris(Time, model3D.Location);
zenith = 90-ApparentSunEl;
azimuth = SunAz;
[x,y,~] = sph2cart((90-azimuth)*deg2radf,(90-zenith)*deg2radf,obj.cbh/cos(zenith*deg2radf));
obj.shadowDisp_d = [-x,-y];
        
frame2 = frame;
for j = 1500/60%10:0.5:50
   r_cloud = 60*j*tan(zenith_cloud*deg2radf);
   sigma = (270-azimuth_cloud)*deg2radf;
   y_dist = r_cloud*sin(sigma);
   x_dist = r_cloud*cos(sigma);
   sensor_location1 = [202.9155,-93.1778]-[x_dist,y_dist];
   sensor_location1_on_image = GetReal2Image(sensor_location1,ocam_model,ecam_model,60*j);
   frame2 = cv.circle(frame2,sensor_location1_on_image,10,'Color',[255,0,0],'Thickness',-1);
end
figure, imshow(frame2);
for i=1:size(bboxes,1)
    % Assume we can approximate the shape of cloud using polytopes
    contourPts = vertcat(contours{1,i}{1,:});
    contourPts = GetImage2Real(contourPts,ocam_model,ecam_model,cbh);
    %x_cloud = [bboxesTemp1(i,1),bboxesTemp1(i,1),bboxesTemp1(i,1)+bboxesTempDim(i,1),bboxesTemp1(i,1)+bboxesTempDim(i,1)];
    %y_cloud = [bboxesTemp1(i,2),bboxesTemp1(i,2)+bboxesTempDim(i,2),bboxesTemp1(i,2)+bboxesTempDim(i,2),bboxesTemp1(i,2)];
    %z_cloud = [cbh, cbh, cbh, cbh]
    x_cloud = contourPts(:,1)';
    y_cloud = contourPts(:,2)';
    z_cloud = cbh*ones(size(contourPts(:,1)'));
    
    % number of points conforming the polynomial
    % approximation of cloud edges
    nclp = numel(z_cloud);
    
    % Plot cloud on 3D map
    cloud = patch ( x_cloud, y_cloud, z_cloud, 0.2 );
    alpha (cloud,0.9);
    
    % Since we know the position, we know x & y coordinates
    % These happen to be the same when projected onto the x-y plane
    % Consider the cloud points determined and time ---- hrs
    time_idx = 1; % 9:00
    %time_idx = 10; % 6:00 pm
    %time_idx = 3 ; % 11:00
    %time_idx = 11; % 7:00 pm
    
    zenith_cloud = zenith(time_idx);
    azimuth_cloud = SunAz(time_idx)';
    
    [x,y,z] = sph2cart((90-azimuth_cloud)*deg2radf,(90-zenith_cloud)*deg2radf,cbh/cos(zenith_cloud*deg2radf));  
    
    % Get the radius of the sun incidence to ground (assuming known CBH)
    %r_cloud = sqrt(x_cloud.^2 + y_cloud.^2);
    r_cloud = cbh*tan(zenith_cloud*deg2radf);
    sigma = (270-azimuth_cloud)*deg2radf;
    y_dist = r_cloud*sin(sigma)*ones(size(y_cloud));
    x_dist = r_cloud*cos(sigma)*ones(size(x_cloud));
     
    shade_orig_x =  x_cloud + x_dist;
    shade_orig_y =  y_cloud + y_dist;
    %  Draw a sun incidence vector projection on ground.
    %hold on,  line ( [shade_orig_x;  x_cloud], ...
    %    [shade_orig_y;  y_cloud], ...
    %    [ones(1,nclp); z_cloud], 'LineStyle','--','Color', 'g', 'LineWidth', 2 )
    
    % % Plot polytope conforming the shadow
    shade = patch ( shade_orig_x, shade_orig_y, 10*ones(1,nclp), 'g' );
end



%%
% Plot the Projected shadow and the sun incidence vector
figure(1)
patch (shade_orig_x, -shade_orig_y,'g') % Note the coordinates for y are inversed here!

% Plot coordinate line on 2D map
hold on, plot([0 0],[-y_lb -y_ub],'m')
hold on, plot([x_ub x_lb],[0 0],'m')

xdata = get(s,'XData');
ydata = get(s,'YData');

xlabel ( '--X axis--' );
ylabel ( '--Y axis--' );
title('2D Shadow Map for CHCRC')

text(5,-200,'North','HorizontalAlignment','left', 'BackgroundColor','w')
text(5,180,'South','HorizontalAlignment','left', 'BackgroundColor','w')
text(-375,10,'West','HorizontalAlignment','left', 'BackgroundColor','w')
text(355,10,'East','HorizontalAlignment','left', 'BackgroundColor','w')

%  Draw a line of the projection of the sun on ground.
%here the z-component is zero
hold on,  line ( [0;  x_sun(time_idx)], ...
    [0;  -y_sun(time_idx)], ...
    'LineStyle','--','Color', 'r', 'LineWidth', 1 )

legend({'Cloud shadow','Incidence vector'}, 'Location','SouthEast')


