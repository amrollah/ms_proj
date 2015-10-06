%% Relative position calculation given two (lat,lon) coordinates.
% phi : lat
% lamb : lon
% var R = 6371; % km
% var phi1 = lat1.toRadians();
% var phi2 = lat2.toRadians();
% var delta_phi = (lat2-lat1).toRadians();
% var delta_lamb = (lon2-lon1).toRadians();

% var a = Math.sin(delta_phi/2) * Math.sin(delta_phi/2) +
%         Math.cos(phi1) * Math.cos(phi2) *
%         Math.sin(delta_lamb/2) * Math.sin(delta_lamb/2);
% var c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1-a));

% var d = R * c;
% var y = Math.sin(lamb2-lamb1) * Math.cos(phi2);
% var x = Math.cos(phi`)*Math.sin(phi2) -
%        Math.sin(phi1)*Math.cos(phi2)*Math.cos(lamb2-lamb1);
% var brng = Math.atan2(y, x).toDegrees();

function [x, y] = GetRelativePosition(loc1,loc2)

%Location.latitude =  47+27/60+33.92/60/60;    %   47°27'33.92" N  Local latitude of CHCRC.C1 roof  (positive if north of equator)
%Location.longitude = 8+16/60+38.50/60/60;     %   8°16'38.50" E   Local Longitude of CHCRC.C1 roof (negative if west)
%Location.altitude =  455;                     %   [m]
%Location.distance = 0;                        %   [m] Distance to the camera
%Location.bearing = 0;                         %   0° 0' 0" Bearing to the camera wrt the North

R = 6371*1000; %km

lat1_r = loc1.latitude*pi/180;
lon1_r = loc1.longitude*pi/180;
lat1 = loc1.latitude;
lon1 = loc1.longitude;

lat2_r = loc2.latitude*pi/180;
lon2_r = loc2.longitude*pi/180;
lat2 = loc2.latitude;
lon2 = loc2.longitude;

delta_lat_r = (lat2-lat1)*pi/180;
delta_lon_r = (lon2-lon1)*pi/180;

a = sin(delta_lat_r/2) * sin(delta_lat_r/2) + ...
    cos(lat1_r) * cos(lat2_r) * ...
    sin(delta_lon_r/2) * sin(delta_lon_r/2);
c = 2 * atan2(sqrt(a),sqrt(1-a));
d = R * c;
y_b = sin(delta_lon_r) * cos(lat2_r);
x_b = cos(lat1_r) * sin(lat2_r) - ...
      sin(lat1_r) * cos(lat2_r) * cos(delta_lon_r);
brng = mod(360 + atan2(y_b,x_b)*180/pi, 360) * pi/180;

x = d * sin(brng);
y = d * cos(brng);
end


