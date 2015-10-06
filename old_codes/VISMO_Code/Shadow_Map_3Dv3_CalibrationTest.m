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

%Conversion to real coordinates using scaled measurements
%L1 = 28.7/1.9*50;  % m
%L2 = 15.3/1.9*50;  % m
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

% figure(1)
% %imagesc(CHCRC,'XData',[x_lb x_ub ],'YData',[y_lb y_ub  ]);
% imagesc(CHCRC,'XData',[x_lb x_ub ],'YData',[-y_ub  -y_lb ]);
% 
% figure(2)
% s = surface('XData',[x_lb x_ub; x_lb x_ub],...
%     'YData',[y_lb y_lb; y_ub y_ub],...
%     'ZData',[0 0; 0 0],'CData',flipdim(CHCRC,1),...
%     'FaceColor','texturemap','EdgeColor','none');
% 
% %text(5,-200,10,'South','HorizontalAlignment','left', 'BackgroundColor','w')
% %text(5,180,10,'North','HorizontalAlignment','left', 'BackgroundColor','w')
% %text(-375,10,10,'West','HorizontalAlignment','left', 'BackgroundColor','w')
% %text(355,10,10,'East','HorizontalAlignment','left', 'BackgroundColor','w')
% 
% xlabel ( '--X axis--' );
% ylabel ( '--Y axis--' );
% zlabel ( '--Z axis--' );
% 
% % Plot coordinate line on 3D map
% hold on, plot([0 0],[y_lb y_ub],'m')
% hold on, plot([x_lb x_ub],[0 0],'m')
% 
% [ x, y, z ] = sphere ( 50 );
% 
% scale_factor  = 200; %Cam_yc;
% 
% x2 = x*scale_factor;
% y2 = y*scale_factor;
% z2 = z*scale_factor;
% 
% z2(z2<0)=0; %form hemisphere
% 
% c = 0.20 * ones ( size ( z2 ) );
% 
% x3 = x*scale_factor*3;
% y3 = y*scale_factor*3;
% z3 = z*scale_factor*3;
% z3(z3<0)=0; %form hemisphere
% c3 = 0.40 * ones ( size ( z3 ) );
% 
% x4 = x./z*(scale_factor*3);
% y4 = y./z*(scale_factor*3);
% z4 = z./z*(scale_factor*3);
% z4(z4<0)=0; %form hemisphere
% c4 = 0.80 * ones ( size ( z4 ) );
% %Plot hemispheric view of fish eye
% esph = surf ( x2, y2, z2, c, 'EdgeColor', 'none');
% alpha(esph,0.1)


%esph3 = surf ( x3, y3, z3, c3, 'EdgeColor', 'none');
%alpha(esph3,0.1)

%esph4 = surf ( x4, y4, z4, c4, 'EdgeColor', 'none');
%alpha(esph4,0.1)

%view(-20,49)


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

% DomeAltitude = 200; %m
% 
% %  Draw the focus, center, and pole points in black.
% hold on,  scatter3 ( [0, 0], [0, 0], [0, DomeAltitude], 'filled', 'k', 'SizeData', 100 );
% 
% %  Draw a line through the focus (also center), and noth pole.
% hold on,  line ( [0,0], [0,0], [0, DomeAltitude], 'Color', 'k', 'LineWidth', 2 )
% 
% %  Draw the image points in red.
% hold on, scatter3 ( x_sun', y_sun', z_sun', 'filled', 'r' );
% 
% %  Draw a line indicating the incidence vector towards the focus.
% hold on,  line ( [zeros(1,nDN);  x_sun], ...
%     [zeros(1,nDN);  y_sun], ...
%     [zeros(1,nDN);  z_sun], 'Color', 'r', 'LineWidth', 1 )
% 
% %  Draw a line of the projection of the sun on ground.
% %here the z-component is zero
% hold on,  line ( [zeros(1,nDN);  x_sun], ...
%     [zeros(1,nDN);  y_sun], ...
%     zeros(2,nDN), 'LineStyle','--','Color', 'r', 'LineWidth', 1)
% 
% %  Draw a line of the sun trajectory -here the z-component is zero
% hold on,  line ( x_sun, ...
%     y_sun, ...
%     z_sun, 'LineStyle','-','Color', 'y', 'LineWidth', 2 );
% 
% title ( 'Stereographic Projection of Sun Trajectory' )



% ----------------------- Shadow Detection ----------------------------
load 'calib_modelGood.mat';
load 'extrinsic_calib_model4.mat';
load 'image_mask.mat';
sunPattern = imread('sunPattern.jpg');

% Assume we know the position of the Cloud
cbh = 1500; % provided by the calibration model, cbh is the assumed height
% Here the cbh has been scaled down for visual inspection purpuses
frame = imread('oimage5.jpeg');
frame(~image_mask) = 0;  % we eliminate the parts that are angularly behind the camera observation point
res = cv.matchTemplate(frame,sunPattern); % find the position of the sun in the current frame through template matching with a previously obtained sun template
[x,y]=ind2sub(size(res),find(res==min(min(res)))); % find the pixel sun position on the image
frame = cv.circle(frame,[y+6,x+2],80,'Color',[0,0,0],'Thickness',-1); % eliminate the sun putting a black circle on top of it.
r_b = double(frame(:,:,1))./double(frame(:,:,3)); % R/B thresholding applied
rgb = double(frame(:,:,1))+double(frame(:,:,2))+double(frame(:,:,3)); % to distinguish other objects like buildings from clouds, a test on brightness
mask = r_b > 0.72 & rgb > 250;
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
centroidsProj = centroidsProj(consider);
bboxes = bboxes(consider);
areas = areas(consider);
contours = contours(consider);
contoursProj = contoursProj(consider);

centroids = vertcat(centroids{:});
centroidsProj = vertcat(centroidsProj{:});

centroidsProjTest=[];
count2 = 1;
range2 = -2:0.1:2;
range3 = -2:0.1:2;
for theta2 = range2
    count3 = 1;
    for theta3 = range3
        contoursProjTest = cellfun(@(x) num2cell(GetImage2RealTest(cell2mat(x(:)),ocam_model,ecam_model,cbh,theta2,theta3),2)',contours,'UniformOutput',false);
        centroidsProjTestPts = cellfun(@(x) (min(cell2mat(x(:)),[],1)+max(cell2mat(x(:)),[],1))/2,contoursProjTest,'UniformOutput',false);
        centroidsProjTestDiff = vertcat(centroidsProjTestPts{:})-centroidsProj;
        centroidsProjTest{count2,count3} = sqrt(centroidsProjTestDiff(:,1).^2 + centroidsProjTestDiff(:,2).^2); 
        count3 = count3+1;
    end
    count2 = count2+1;
end

%%
for count2 = 1:size(centroidsProjTest,1)
    for count3 = 1:size(centroidsProjTest,2)
        centroidsProjTestMean(count2,count3) = mean(centroidsProjTest{count2,count3});
        centroidsProjTestMin(count2,count3) = min(centroidsProjTest{count2,count3});
        centroidsProjTestMax(count2,count3) = max(centroidsProjTest{count2,count3});
    end
end
centroidsProjTestMean(centroidsProjTestMean > 10^8) = 10^6;
centroidsProjTestMax(centroidsProjTestMax > 10^8) = 10^6;

figure, surf(range2,range3,centroidsProjTestMin);
xlabel('Roll (degrees)','FontSize',18), ylabel('Pitch (degrees)','FontSize',18), title('Changes in projection with assumed rotations of the camera (meters)','FontSize',18), set(gca,'FontSize',18),axis equal, colorbar;
figure, surf(range2,range3,centroidsProjTestMean);
xlabel('Rotations in degrees around x axis'), ylabel('Rotations in degrees around y axis'), title('Mean changes in cloud position with assumed rotations of the camera'), colorbar;
figure, surf(range2,range3,centroidsProjTestMax);
xlabel('Rotations in degrees around x axis'), ylabel('Rotations in degrees around y axis'), title('Maximum changes in cloud position with assumed rotations of the camera'), colorbar;

%%
% bboxes = vertcat(bboxes{:}); % put all the boxes in a 2D array, so that the annotation algorithm can display them
% bboxesTemp1 = GetImage2Real(bboxes(:,1:2),ocam_model,ecam_model,cbh); % upper corner of bboxes in projected world
% bboxesTemp2 = GetImage2Real(bboxes(:,1:2)+bboxes(:,3:4),ocam_model,ecam_model,cbh); % lower corner of bboxes in projected world
% bboxesTempDim = abs(bboxesTemp2 - bboxesTemp1); % dimensions of bboxes (in a sense)
% for i=1:size(bboxes,1)
%     % Assume we can approximate the shape of cloud using polytopes
%     contourPts = vertcat(contours{1,i}{1,:});
%     contourPts = GetImage2Real(contourPts,ocam_model,ecam_model,cbh);
%     count2 = 1;
%     for theta2 = -10:0.5:10
%         count3 = 1;
%         for theta3 = -10:0.5:10
%             centroidsProjTest{count2,count3} = GetImage2RealTest(centroids(:,1:2),ocam_model,ecam_model,cbh,theta2,theta3)-centroidsProj;
%             count3 = count3+1;
%         end
%         count2 = count2+1;
%     end
%     %x_cloud = [bboxesTemp1(i,1),bboxesTemp1(i,1),bboxesTemp1(i,1)+bboxesTempDim(i,1),bboxesTemp1(i,1)+bboxesTempDim(i,1)];
%     %y_cloud = [bboxesTemp1(i,2),bboxesTemp1(i,2)+bboxesTempDim(i,2),bboxesTemp1(i,2)+bboxesTempDim(i,2),bboxesTemp1(i,2)];
%     %z_cloud = [cbh, cbh, cbh, cbh]
%     x_cloud = contourPts(:,1)';
%     y_cloud = contourPts(:,2)';
%     z_cloud = cbh*ones(size(contourPts(:,1)'));
%     
%     % number of points conforming the polynomial
%     % approximation of cloud edges
%     nclp = numel(z_cloud);
%     
%     % Plot cloud on 3D map
%     cloud = patch ( x_cloud, y_cloud, z_cloud, 0.2 );
%     alpha (cloud,0.9);
%     
%     % Since we know the position, we know x & y coordinates
%     % These happen to be the same when projected onto the x-y plane
%     % Consider the cloud points determined and time ---- hrs
%     time_idx = 1; % 9:00
%     %time_idx = 10; % 6:00 pm
%     %time_idx = 3 ; % 11:00
%     %time_idx = 11; % 7:00 pm
%     
%     zenith_cloud = zenith(time_idx);
%     azimuth_cloud = SunAz(time_idx)';
%     
%     [x,y,z] = sph2cart((90-azimuth_cloud)*deg2radf,(90-zenith_cloud)*deg2radf,cbh/cos(zenith_cloud*deg2radf));
%     
%     
%     % Get the radius of the sun incidence to ground (assuming known CBH)
%     %r_cloud = sqrt(x_cloud.^2 + y_cloud.^2);
%     r_cloud = cbh*tan(zenith_cloud*deg2radf);
%     sigma = (270-azimuth_cloud)*deg2radf;
%     y_dist = r_cloud*sin(sigma)*ones(size(y_cloud));
%     x_dist = r_cloud*cos(sigma)*ones(size(x_cloud));
%     
%     shade_orig_x =  x_cloud + x_dist;
%     shade_orig_y =  y_cloud + y_dist;
%     %  Draw a sun incidence vector projection on ground.
%     %hold on,  line ( [shade_orig_x;  x_cloud], ...
%     %    [shade_orig_y;  y_cloud], ...
%     %    [ones(1,nclp); z_cloud], 'LineStyle','--','Color', 'g', 'LineWidth', 2 )
%     
%     % % Plot polytope conforming the shadow
%     shade = patch ( shade_orig_x, shade_orig_y, zeros(1,nclp), 'g' );
% end
% 
% 
% 
% %%
% % Plot the Projected shadow and the sun incidence vector
% figure(1)
% patch (shade_orig_x, -shade_orig_y,'g') % Note the coordinates for y are inversed here!
% 
% % Plot coordinate line on 2D map
% hold on, plot([0 0],[-y_lb -y_ub],'m')
% hold on, plot([x_ub x_lb],[0 0],'m')
% 
% xdata = get(s,'XData');
% ydata = get(s,'YData');
% 
% xlabel ( '--X axis--' );
% ylabel ( '--Y axis--' );
% title('2D Shadow Map for CHCRC')
% 
% text(5,-200,'North','HorizontalAlignment','left', 'BackgroundColor','w')
% text(5,180,'South','HorizontalAlignment','left', 'BackgroundColor','w')
% text(-375,10,'West','HorizontalAlignment','left', 'BackgroundColor','w')
% text(355,10,'East','HorizontalAlignment','left', 'BackgroundColor','w')
% 
% %  Draw a line of the projection of the sun on ground.
% %here the z-component is zero
% hold on,  line ( [0;  x_sun(time_idx)], ...
%     [0;  -y_sun(time_idx)], ...
%     'LineStyle','--','Color', 'r', 'LineWidth', 1 )
% 
% legend({'Cloud shadow','Incidence vector'}, 'Location','SouthEast')
% 
% 
% 
% 
% % %-----------------------------------
% % % Load VISMO image
% %
% % %load picture
% % %figure(3), imshow(H);
% % figure, vismopic = imagesc(H), axis equal,
% % %rotate(vismopic, [1 0],20)
% %
% % [x_pix_length, y_pix_length, rgb] = size(H)
% %
% % % Assume the focus is center of the image
% % x_pix_c = x_pix_length/2; % Center
% % y_pix_c = y_pix_length/2; % Center
% %
% % hold on, plot(x_pix_c,y_pix_c, 'LineWidth',4,'Color','g')
% %
% 
% 
% 
% 
% 
% 
% 
% %addpath('X:\ABB Projects\Solar\Cloud Tracking\mexopencv-master')
% %rotated_img=cv.rotateImage(H,10);
% %help cv.imread
% 
% 
% % rotatedImg = imrotate(H,30); % Try varying the angle, theta.
% % figure, imshow(rotatedImg)
% % figure, imagesc(rotatedImg), axis equal,
% 

%%
%47°27'36.0"N 8°17'21.6"E
%47°27'35.8"N 8°16'37.5"E

% ref 47.459954, 8.277085
% 15 km East 47.459201, 8.458732
% 15 km South 47.336961, 8.279518
% 15 km West 47.460826, 8.094124
% 15 km North 47.583010, 8.276428
% 5 km North 47.501663, 8.277458
% 5 km East 47.460130, 8.337539
% 5 km West 47.461987, 8.215317
% 5 km South 47.419028, 8.278145



% alternative locations 
% loc E : 47°27'37.1"N 8°23'37.4"E
% loc N : 47°32'01.0"N 8°16'15.0"E
% loc W : 47°28'08.6"N 8°07'27.9"E
% loc S : 47°21'44.8"N 8°17'39.4"E
% loc C : 47°27'38.0"N 8°16'36.5"E
% ln.latitude =  47+32/60+1.0/60/60;    %  47°32'01.0"N   Local latitude 
% ln.longitude = 8+16/60+15.0/60/60;     %   8°16'15.0"E   Local Longitude
% ln.altitude =  0;
% ls.latitude =  47+21/60+44.8/60/60;    %  47°21'44.8"N   Local latitude 
% ls.longitude = 8+17/60+39.4/60/60;     %   8°17'39.4"E   Local Longitude 
% ls.altitude =  0;
% le.latitude =  47+27/60+37.1/60/60;    %  47°27'37.1"N   Local latitude 
% le.longitude = 8+23/60+37.4/60/60;     %   8°23'37.4"E   Local Longitude
% le.altitude =  0;
% lw.latitude =  47+28/60+8.6/60/60;    %  47°28'08.6"N   Local latitude 
% lw.longitude = 8+7/60+27.9/60/60;     %   8°07'27.9"E   Local Longitude 
% lw.altitude =  0;

% 15 km around 47.459954, 8.277085
dist = 3e4;
ln.latitude =  47.583010;    
ln.longitude = 8.276428;     
ln.altitude =  0;
ls.latitude =  47.336961;    
ls.longitude = 8.279518;     
ls.altitude =  0;
le.latitude =  47.459201;    
le.longitude = 8.458732;     
le.altitude =  0;
lw.latitude =  47.460826;    
lw.longitude = 8.094124; 
lw.altitude =  0;

% 5 km around 47.459954, 8.277085
% dist = 1e4;
% ln.latitude =  47.501663;    
% ln.longitude = 8.277458;    
% ln.altitude =  0;
% ls.latitude =  47.419028;   
% ls.longitude = 8.278145;   
% ls.altitude =  0;
% le.latitude =  47.460130;   
% le.longitude = 8.337539;     
% le.altitude =  0;
% lw.latitude =  47.461987;    
% lw.longitude = 8.215317;     
% lw.altitude =  0;

cbh_def  = 1500;
lref.latitude =  47.459954;    
lref.longitude = 8.277085;     
lref.altitude =  0;


lnw.latitude = ln.latitude;
lnw.longitude = lw.longitude;
lnw.altitude = 0;
lse.latitude = ls.latitude;
lse.longitude = le.longitude;
lse.altitude = 0;


start_date = datenum(2013,07,12);
dt = 12/24;
DN = dt +  start_date;
nDN= numel(DN);

UTC = 2; %
Time = pvl_maketimestruct(DN, UTC);
% Pre-calculations & conversion factors
rad2degf = 180/pi;
deg2radf = pi/180;

d_lat = (lnw.latitude - lse.latitude)/100;
d_long = (lse.longitude - lnw.longitude)/100;
d_dist = dist/100;
% form a grid of locations:

[SunAz, ~, ApparentSunEl, ~]=pvl_ephemeris(Time, lref);
zenithref = 90-ApparentSunEl;
azimuthref = SunAz;
[xref,yref,~] = sph2cart((90-azimuthref)*deg2radf,(90-zenithref)*deg2radf,cbh_def/cos(zenithref*deg2radf));

for i = 1:100
    for j = 1:100
        temp.latitude = lse.latitude + (i-1)*d_lat;
        temp.longitude = lnw.longitude + (j-1)*d_long;
        temp.altitude = 0;
        loc_grid(i,j) = temp;
        
        [SunAz, ~, ApparentSunEl, ~]=pvl_ephemeris(Time, temp);
        zenith = 90-ApparentSunEl;
        azimuth = SunAz;
        [x,y,~] = sph2cart((90-azimuth)*deg2radf,(90-zenith)*deg2radf,cbh_def/cos(zenith*deg2radf));
        
        sun_dir(i,j).x = x;
        sun_dir(i,j).y = y;
        
        deg_gridy(i,j) = temp.latitude - lref.latitude;
        deg_gridx(i,j) = temp.longitude - lref.longitude;
        loc_gridy(i,j) = -dist/2  + (i-1)*d_dist;
        loc_gridx(i,j) = -dist/2  + (j-1)*d_dist;
        sun_dirx(i,j) = x-xref;
        sun_diry(i,j) = y-yref;
        sun_zen(i,j) = zenith-zenithref;
        sun_az(i,j) = azimuth-azimuthref;
    end
end

%%
figure, scatter(sun_dirx(:),sun_diry(:));
figure, scatter(sun_zen(:),sun_az(:));


figure, surf(loc_gridx,loc_gridy,abs(sun_zen));
figure, surf(loc_gridx,loc_gridy,abs(sun_az));
figure, surf(loc_gridx,loc_gridy,abs(sun_dirx));
figure, surf(loc_gridx,loc_gridy,abs(sun_diry));
figure, surf(loc_gridx,loc_gridy,sqrt(sun_dirx.^2+sun_diry.^2));
figure, scatter(sqrt(loc_gridx(:).^2+loc_gridy(:).^2),sqrt(sun_dirx(:).^2+sun_diry(:).^2));

figure, surf(loc_gridx(45:55,45:55),loc_gridy(45:55,45:55),sun_zen(45:55,45:55));
figure, surf(loc_gridx(45:55,45:55),loc_gridy(45:55,45:55),sun_az(45:55,45:55));
figure, surf(loc_gridx(25:75,25:75),loc_gridy(25:75,25:75),sun_dirx(25:75,25:75));
figure, surf(loc_gridx(25:75,25:75),loc_gridy(25:75,25:75),sun_diry(25:75,25:75));

%% sample
loc1.latitude = 47+27/60+35.8/60/60;
loc1.longitude = 8+16/60+37.5/60/60;
loc1.altitude = 0;

[SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(Time, loc1);
zenithl1 = 90-ApparentSunEl;
azimuthl1 = SunAz;
[xl1,yl1,~] = sph2cart((90-azimuthl1)*deg2radf,(90-zenithl1)*deg2radf,1/cos(zenithl1*deg2radf));

loc2.latitude = 47+27/60+36.0/60/60;
loc2.longitude = 8+17/60+21.6/60/60;
loc2.altitude = 0;

[SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(Time, loc2);
zenithl2 = 90-ApparentSunEl;
azimuthl2 = SunAz;
[xl2,yl2,~] = sph2cart((90-azimuthl2)*deg2radf,(90-zenithl2)*deg2radf,1/cos(zenithl2*deg2radf));