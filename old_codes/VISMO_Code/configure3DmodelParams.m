%% generic parameters of model3d
Location.latitude =  47.459223;    %   47°27'33.92" N  Local latitude of CHCRC.C1 roof  (positive if north of equator) 47.459223 47+27/60+33.92/60/60
Location.longitude = 8.277023;     %   8°16'38.50" E   Local Longitude of CHCRC.C1 roof (negative if west) 8.277023 8+16/60+38.50/60/60
Location.altitude =  455;                     %   [m]
Location.distance = 0;                        %   [m] Distance to the camera
Location.bearing = 0;                         %   0° 0' 0" Bearing to the camera wrt the North

%Sensor on the roof top.
%sensor location 2 47.458385, 8.279722
SensorLocation.latitude =  47.459519;    %   47° 27' 34.27" N  Local latitude of CHCRC.C1 roof  (positive if north of equator) 47.459519 47+27/60+34.27/60/60
SensorLocation.longitude = 8.277495;     %   8° 16' 38.98" E    Local Longitude of CHCRC.C1 roof (negative if west) 8.277495 8+16/60+38.98/60/60
SensorLocation.altitude =  455;                     %   [m]
SensorLocation.distance = 48;                       %   [m] Distance to the camera
SensorLocation.bearing = 47+9/60+10/60/60;          %   47° 09' 10" Bearing to the camera wrt the North 47+9/60+10/60/60

%Sensor at the parking place.
SensorLocation2.latitude = 47.458385; 
SensorLocation2.longitude = 8.279722;
SensorLocation2.altitude = 455;
SensorLocation2.distance = -1;
SensorLocation2.bearing = -1;

% the other sensors are the three sensors surrounding the Baden-Segelhof
% office in the map. You can find the exact correspondance by looking at
% the positions
% http://www.wunderground.com/wundermap/?lat=47.47999954&lon=8.39999962&zoom=8&pin=Laegern%2c%20&rad=0&rad.type=00Q&wxsn=1&svr=0&cams=0&sat=0&riv=0&mm=0&hur=0

%47.481 8.293 http://www.wunderground.com/personal-weather-station/dashboard?ID=IAARGAUB5#history/s20140612/e20140612/mdaily
SensorLocation3.latitude = 47.481;
SensorLocation3.longitude = 8.293;
SensorLocation3.altitude = 455;
SensorLocation3.distance = -1;
SensorLocation3.bearing = -1;

%47.411 8.316 http://www.wunderground.com/personal-weather-station/dashboard?ID=IAARGAUU2#history/s20130812/e20130812/mdaily
SensorLocation4.latitude = 47.411;
SensorLocation4.longitude = 8.316;
SensorLocation4.altitude = 455;
SensorLocation4.distance = -1;
SensorLocation4.bearing = -1;

%47.505 8.248 http://www.wunderground.com/personal-weather-station/dashboard?ID=IAARGAUU2#history/s20130812/e20130812/mdaily
SensorLocation5.latitude = 47.505;
SensorLocation5.longitude = 8.248;
SensorLocation5.altitude = 455;
SensorLocation5.distance = -1;
SensorLocation5.bearing = -1;

model3D.Location = Location;
model3D.SensorLocation = SensorLocation;
model3D.SensorLocation2 = SensorLocation2;
model3D.SensorLocation3 = SensorLocation3;
model3D.SensorLocation4 = SensorLocation4;
model3D.SensorLocation5 = SensorLocation5;
model3D.UTC = 2;

%model3D.Sen_xs = SensorLocation.distance*sin(SensorLocation.bearing*pi/180);
%model3D.Sen_ys = SensorLocation.distance*cos(SensorLocation.bearing*pi/180);
[model3D.Sen_xs, model3D.Sen_ys] = GetRelativePosition(Location,SensorLocation);
[model3D.Sen2_xs, model3D.Sen2_ys] = GetRelativePosition(Location,SensorLocation2);
[model3D.Sen3_xs, model3D.Sen3_ys] = GetRelativePosition(Location,SensorLocation3);
[model3D.Sen4_xs, model3D.Sen4_ys] = GetRelativePosition(Location,SensorLocation4);
[model3D.Sen5_xs, model3D.Sen5_ys] = GetRelativePosition(Location,SensorLocation5);

% ************** the way to do the occlusion and irradiance predictions for
% a set of points is as follows: 
% - a location should be added just like the ones above.
% - the position of the added location should be calculated wrt the camera
% as done just above.
% - the x and y coordinates should be added to the below list.
% - inside the code at the part marked by the ***** the below lists should
% be used to replace model3D.Sen_xs.
model3D.all_xs = [model3D.Sen_xs; model3D.Sen2_xs; model3D.Sen3_xs; model3D.Sen4_xs; model3D.Sen5_xs];
model3D.all_ys = [model3D.Sen_ys; model3D.Sen2_ys; model3D.Sen3_ys; model3D.Sen4_ys; model3D.Sen5_ys];
%% for the image 'CHCRC_Bigger.jpg'
L1 = 40.9/2*500;  % m
L2 = 21.1/2*500;  % m

% Camera Coordinates on map
model3D.Cam_xc = 20.5/2*500;  %m
model3D.Cam_yc = 10.1/2*500;  %m

% Remaping cordinates,  camera as origin
model3D.x_lb =    - model3D.Cam_xc;
model3D.x_ub = L1 - model3D.Cam_xc;

model3D.y_lb =    - model3D.Cam_yc;
model3D.y_ub = L2 - model3D.Cam_yc;

model3D.all = [model3D.x_lb,model3D.x_ub,model3D.y_lb,model3D.y_ub,model3D.Cam_xc,model3D.Cam_yc,model3D.Sen_xs,model3D.Sen_ys];
model3D.fileName = 'CHCRC_Bigger.jpg';

cd Utils/Data
save('model3D_Bigger', 'model3D');
cd ../../
%% for the image 'CHCRC_Bigger2.jpg'
L1 = 40.9/2*1000;  % m
L2 = 21.1/2*1000;  % m

% Camera Coordinates on map
model3D.Cam_xc = 20.5/2*1000;  %m
model3D.Cam_yc = 10.1/2*1000;  %m

% Remaping cordinates,  camera as origin
model3D.x_lb =    - model3D.Cam_xc;
model3D.x_ub = L1 - model3D.Cam_xc;

model3D.y_lb =    - model3D.Cam_yc;
model3D.y_ub = L2 - model3D.Cam_yc;

model3D.all = [model3D.x_lb,model3D.x_ub,model3D.y_lb,model3D.y_ub,model3D.Cam_xc,model3D.Cam_yc,model3D.Sen_xs,model3D.Sen_ys];
model3D.fileName = 'CHCRC_Bigger2.jpg';

cd Utils/Data
save('model3D_Bigger2', 'model3D');
cd ../../
%% for the image 'CHCRC_Bigger3.jpg'
L1 = 40.9/2.5*5000;  % m
L2 = 21.1/2.5*5000;  % m

% Camera Coordinates on map
model3D.Cam_xc = 20.5/2.5*5000;  %m
model3D.Cam_yc = 10.1/2.5*5000;  %m

% Remaping cordinates,  camera as origin
model3D.x_lb =    - model3D.Cam_xc;
model3D.x_ub = L1 - model3D.Cam_xc;

model3D.y_lb =    - model3D.Cam_yc;
model3D.y_ub = L2 - model3D.Cam_yc;

model3D.all = [model3D.x_lb,model3D.x_ub,model3D.y_lb,model3D.y_ub,model3D.Cam_xc,model3D.Cam_yc,model3D.Sen_xs,model3D.Sen_ys];
model3D.fileName = 'CHCRC_Bigger3.jpg';

cd Utils/Data
save('model3D_Bigger3', 'model3D');
cd ../../

%% for the image 'KACST.bmp'
L1 = 28.8/1.8*500;  % m
L2 = 15.2/1.8*500;  % m

% Camera Coordinates on map
model3D.Cam_xc = 13.1/1.8*500;  %m
model3D.Cam_yc = 5.5/1.8*500;  %m

% Remaping cordinates,  camera as origin
model3D.x_lb =    - model3D.Cam_xc;
model3D.x_ub = L1 - model3D.Cam_xc;

model3D.y_lb =    - model3D.Cam_yc;
model3D.y_ub = L2 - model3D.Cam_yc;

model3D.all = [model3D.x_lb,model3D.x_ub,model3D.y_lb,model3D.y_ub,model3D.Cam_xc,model3D.Cam_yc,model3D.Sen_xs,model3D.Sen_ys];
model3D.fileName = 'KACST.bmp';

cd Utils/Data
save('model3D_KACST', 'model3D');
cd ../../
