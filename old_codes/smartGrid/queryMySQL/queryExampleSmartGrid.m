% This file identifes the ID_links of the available sensors 

% 21.06.2013
% Luis Dominguez
% Control & Optimization
% ABB Schweiz
%-------------------------------------------------------------------------


clear all
clf

%% set path
addpath(fullfile(pwd, 'src'));
javaaddpath('lib/mysql-connector-java-5.1.6/mysql-connector-java-5.1.6-bin.jar');

%% import classes
import edu.stanford.covert.db.MySQLDatabase;

%% create database connection
% Smart meters:
% http://10.41.94.130   - Works
% http://10.41.94.131   - Doesn't Work

%db = MySQLDatabase('covertlab.stanford.edu', 'test', 'test', 'test');
db = MySQLDatabase('10.41.94.39', 'smartgrid', 'root', 'smartgrid');




%% send, retrieve blob using sql

% ----------------------------Available Sensors of Wheather Station ------
% 0.3;    weather station;
% 0.3.1;  brightness right;
% 0.3.2;  brightness center;
% 0.3.3;  brightness left;
% 0.3.4;  outside temperature;
% 0.3.5;  wind speed;
% 0.3.11; brightness right threshold;
% 0.3.12; center brightness threshold;
% 0.3.13; brightness threshold left;
% 0.3.15; wind speed threshold;
% 0.3.16; rain threshold;
% 0.3.21, change the brightness threshold lower limit;
% 0.3.22, change the brightness threshold upper limit;
% 0.3.51; facade A Nord-Ost approach position;
% 0.3.52; facade A Nord-Ost approach lamella;
% 0.3.53; facade A Nord-Ost sun;
% 0.3.61; approach facade B South East position;
% 0.3.62; facade B Move South East lamella;
% 0.3.63; facade B South East sun;
% 0.3.71; facade C South-West approach position;
% 0.3.72; facade C South-West Approach lamella;
% 0.3.73; facade C South-western sun;
% 0.3.81; facade D North West approach position;
% 0.3.82; facade D North West Approach lamella;
% 0.3.83; facade D North West sun;
% 0.3.100; Time;
% 0.3.101; date;
% 0.3.102; time requirement;
% 0.3.103; status WZ/S 1.1;
% 0.3.104; sensor failure WZ/S 1.1;
% 0.3.105; no time synchronization;
% 0.3.200; pyranometer 0-2000W/m2;
% 0.3.201; WS/S Status byte;

% To get the ID_link from Pyranometer type:
%db.prepareStatement('select * from id_link where UniqueID="0.3.200"');
%id_Pyr = db.query();

% To get the ID_link from Wind Speed type:
db.prepareStatement('select * from id_link where UniqueID="0.3.5"');
id_WindSpeed = db.query();

% To get the ID_link from Outside Temerature type:
db.prepareStatement('select * from id_link where UniqueID="0.3.4"');
id_AmbTemp = db.query();

% To get the ID_link from Wind speed threshold type:
db.prepareStatement('select * from id_link where UniqueID="0.3.15"');
id_WindSpeedThreshold = db.query();

% To get the ID_link from Weather Sensor Status type:
db.prepareStatement('select * from id_link where UniqueID="0.3.201"');
id_WeatherSensorStatus = db.query();

% To get the ID_link from status WZ/S (weather Unit) Status type:
db.prepareStatement('select * from id_link where UniqueID="0.3.103"');
id_WeatherUnitStatus = db.query();

% To get the ID_link from Brightness Left type:
db.prepareStatement('select * from id_link where UniqueID="0.3.3"');
id_BirghtLeft = db.query();

% To get the ID_link from Brightness Right type:
db.prepareStatement('select * from id_link where UniqueID="0.3.1"');
id_BirghtRight = db.query();

% To get the ID_link from Brightness Center type:
db.prepareStatement('select * from id_link where UniqueID="0.3.2"');
id_BirghtCenter = db.query();

% select * from id_link where uniqueID like "ABB Electricity Meter 00411665.%";

% +-------+-----------------------------------------------------+
% | id    | uniqueID                                            |
% +-------+-----------------------------------------------------+
% | 16532 | ABB Electricity Meter 00411665.Current L1           |
% | 16533 | ABB Electricity Meter 00411665.Current L2           |
% | 16534 | ABB Electricity Meter 00411665.Current L3           |
% | 16524 | ABB Electricity Meter 00411665.Energy Active export |
% | 16523 | ABB Electricity Meter 00411665.Energy Active import |
% | 16535 | ABB Electricity Meter 00411665.Frequency            |
% | 16525 | ABB Electricity Meter 00411665.Power Active         |
% | 16526 | ABB Electricity Meter 00411665.Power Active L1      |
% | 16527 | ABB Electricity Meter 00411665.Power Active L2      |
% | 16528 | ABB Electricity Meter 00411665.Power Active L3      |
% | 16536 | ABB Electricity Meter 00411665.Power Factor         |
% | 16529 | ABB Electricity Meter 00411665.Voltage L1           |
% | 16530 | ABB Electricity Meter 00411665.Voltage L2           |
% | 16531 | ABB Electricity Meter 00411665.Voltage L3           |
% +-------+-----------------------------------------------------+


%select * from id_16663 where timestamp between '2013-05-01 05:00' and '2013-05-13 21:00'  limit 1000000;
%select * from id_16663 where timestamp between '2013-05-15 00:00' and '2013-05-15 23:59'  limit 1000000;

meas = 10 % Number of measurments required
lim = num2str (meas);

% ---------------------  Atmospheric Readings ----------------------------
%Weather Sensor Status
db.prepareStatement(['select * from id_16664 order by timestamp desc limit ',lim]);
WeatherSensorSatus_result = db.query();

WeatherSensorSatus_result.timestamp(1);
WeatherSensorSatus_result.value(1);


%Weather  Status
db.prepareStatement(['select * from id_16660 order by timestamp desc limit ',lim]);
WeatherUnitStatus_result = db.query();

WeatherUnitStatus_result.timestamp(1);
WeatherUnitStatus_result.value(1);


% Pyranometer
db.prepareStatement(['select * from id_16663 order by timestamp desc limit ',lim]);
%db.prepareStatement('select * from id_16663 where timestamp between "2013-06-10 00:00" and "2013-06-10 23:59"  limit 1000000');

GHI_result = db.query();
GHI_result.timestamp(1);
GHI_result.value(1) * 100; %Multiply by 100 to match KNX display 

%from past to now
dt_GHI = datenum(GHI_result.timestamp(end:-1:1));
GHI = GHI_result.value (end:-1:1) * 100;  %Multiply by 100 to match KNX display 

datestr(dt_GHI)
GHI

return

% % %Wind Speed
db.prepareStatement(['select * from id_16637 order by timestamp desc limit ',lim]);
WindSpeed_result = db.query();

WindSpeed_result.timestamp(1);  % Wind does not work !!
WindSpeed_result.value(1);

% % Wind Speed Threshold
% db.prepareStatement('select * from id_16641 order by timestamp desc limit ',lim]); 
% WindSpeedThreshold_result = db.query();
%  
% WindSpeedThreshold_result.timestamp(1)  
% WindSpeedThreshold_result.value(1)

%Temperature
db.prepareStatement(['select * from id_16636 order by timestamp desc limit ',lim]); 
Temp_result = db.query();

Temp_result.timestamp(1);  % Temp & Pyr are almost syncronized - Thy have the same dates
Temp_result.value(1);

%from past to now
dt_Temp = datenum(Temp_result.timestamp(end:-1:1));
Temp = Temp_result.value (end:-1:1);  %


% Brightness Left
 db.prepareStatement(['select * from id_16635 order by timestamp desc limit ',lim]); 
BrightLeft_result = db.query();

BrightLeft_result.timestamp(1); 
BrightLeft_result.value(1) /1000;  % divide by 1000 to get Lux units

% Brightness Center
db.prepareStatement(['select * from id_16634 order by timestamp desc limit ',lim]); 
BrightCenter_result = db.query();

BrightCenter_result.timestamp(1) 
BrightCenter_result.value(1) /1000;  % divide by 1000 to get Lux units

% Brightness Right
db.prepareStatement(['select * from id_16633 order by timestamp desc limit ',lim]); 
BrightRight_result = db.query();

BrightRight_result.timestamp(1); 
BrightRight_result.value(1) /1000;  % divide by 1000 to get Lux units


% ----------------------  Power Generation Readings -----------------------
% Power
db.prepareStatement(['select * from id_16525 order by timestamp desc limit ',lim]);
PV_result = db.query();

PV_result.timestamp(1);  % is the latests 
PV_result.value(1);
 
%from past to now
dt_PV = datenum(PV_result.timestamp(end:-1:1));
P_PV  = -1*PV_result.value (end:-1:1);  % multiply by -1 to get positive values in Power


% ----------------------- Resample Procedure ----------------------------
% Variables have different time stamps
% Resampling variables / to make them all uniform

now_vec = datevec(now); 
%now_vec = datevec(dt_GHI(end))
%datestr(now_vec)

day_now = datenum(now_vec(1:3));
%datestr(day_now) % debug
hr_now  =  datenum(now_vec(4));
min_now = datenum(now_vec(5)) ;

% 10 min interval 1min resolution
last_ten_min  = 1 /24/60:1 /24/60:10 /24/60 ;

% Up to minuntes
now_time = day_now + hr_now/24 + min_now/24/60; 
%datestr(now_time) % Debug

% Past window
dt_past_window = now_time - last_ten_min';
%datestr(dt_past_window) % Debug

% Resample all variables using the past window
GHI_past  = interp1(dt_GHI,GHI,dt_past_window);
Temp_past = interp1(dt_Temp,Temp,dt_past_window);
P_PV_past = interp1(dt_PV,P_PV,dt_past_window);



%-------------- Plot Figures -------------------------------------------
figure(1)
%GHI
subplot(3,1,1)
plot(dt_GHI,GHI), 
hold on, plot(dt_past_window,GHI_past,'m'), hold off 
title(['Measurement History of GHI for ', datestr(floor(dt_GHI(end)))])
ylabel('W/m^2')
xlim([dt_GHI(1) - 1/24/60, dt_GHI(end) + 1/24/60  ])
datetick('x','HH:MM','keeplimits')

% Temperature
subplot(3,1,2)
plot(dt_Temp,Temp),
ylabel ('°C')
title(['Measurement History of Temperature for ',datestr(floor(dt_GHI(end)))])
xlim([dt_GHI(1) - 1/24/60, dt_GHI(end) + 1/24/60  ])
datetick('x','HH:MM','keeplimits')

%Power
subplot(3,1,3)
plot(dt_PV,P_PV), 
ylabel ('W')
title(['Measurement History of Generated Power for ', datestr(floor(dt_PV(end)))])
xlim([dt_GHI(1) - 1/24/60, dt_GHI(end) + 1/24/60])
datetick('x','HH:MM','keeplimits')
%set(gca,'XTick',0:0.2:1); %Set the XTicks values
%datetick('x','HH','keepticks')
%xlim([0 1])
%ylim([0 1.2e3])






%------------- Past Window of 10 min ------------------------------------
figure(2)
subplot(3,1,1)
xlim([dt_past_window(end) - 1/24/60, dt_past_window(1) + 1/24/60  ])
plot(dt_past_window,GHI_past)
title(['10 min Past Window of GHI for ', datestr(floor(dt_GHI(end)))])
ylabel('W/m^2')
datetick('x','HH:MM','keeplimits')

% Temperature
subplot(3,1,2)
xlim([dt_past_window(end) - 1/24/60, dt_past_window(1) + 1/24/60  ])
plot(dt_past_window,Temp_past),
ylabel ('°C')
title(['10 min Past Window of Temperature for ',datestr(floor(dt_GHI(end)))])
datetick('x','HH:MM','keeplimits')

%Power
subplot(3,1,3)
xlim([dt_past_window(end) - 1/24/60, dt_past_window(1) + 1/24/60  ])
plot(dt_past_window,P_PV_past), 
ylabel ('W')
title(['10 min Past Window of Power for ', datestr(floor(dt_PV(end)))])
datetick('x','HH:MM','keeplimits')
%set(gca,'XTick',0:0.2:1); %Set the XTicks values
%datetick('x','HH','keepticks')
%xlim([0 1])
%ylim([0 1.2e3])






%db.close();  %Remember to close the db when finish !!!!
