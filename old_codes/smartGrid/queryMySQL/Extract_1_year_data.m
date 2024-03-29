% This file queries the the SQL database generated by the KNX agents

% 24.06.2013
% Luis Dominguez
% Control & Optimization
% ABB Schweiz
%-------------------------------------------------------------------------

clear all
clf
clc

%% set path
%addpath(fullfile(pwd, 'src'));
%addpath(fullfile(pwd));
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
meas1 = 6500 ; % Number of measurments required for long history
lim = num2str (meas1);

% ---------------------  Atmospheric Readings ----------------------------
%Weather Sensor Status
% db.prepareStatement(['select * from id_16664 order by timestamp desc limit ',lim]);
% WeatherSensorSatus_result = db.query();

% WeatherSensorSatus_result.timestamp(1);
% WeatherSensorSatus_result.value(1);

%Weather  Status
% db.prepareStatement(['select * from id_16660 order by timestamp desc limit ',lim]);
% WeatherUnitStatus_result = db.query();

% WeatherUnitStatus_result.timestamp(1);
% WeatherUnitStatus_result.value(1);


% Pyranometer
%db.prepareStatement(['select * from id_16663 order by timestamp desc limit ',lim]);
db.prepareStatement('select * from id_16663 where timestamp between "2012-06-01 00:00" and "2013-07-02 23:59"  limit 10000000000000');

GHI_result = db.query();

%from past to now
dt_GHI = datenum(GHI_result.timestamp(end:-1:1) );
GHI = GHI_result.value(end:-1:1)  * 100;  %Multiply by 100 to match KNX display 


% % %Wind Speed
%db.prepareStatement(['select * from id_16637 order by timestamp desc limit ',lim]);
% db.prepareStatement('select * from id_16637 where timestamp between "2012-06-010 00:00" and "2013-07-02 23:59"  limit 1e8');
% WindSpeed_result = db.query(); % Wind does not work !!


% % Wind Speed Threshold
% db.prepareStatement('select * from id_16641 order by timestamp desc limit ',lim]); 
% WindSpeedThreshold_result = db.query();
%  
% WindSpeedThreshold_result.timestamp(1)  
% WindSpeedThreshold_result.value(1)

%Temperature
%db.prepareStatement(['select * from id_16636 order by timestamp desc limit ',lim]); 
db.prepareStatement('select * from id_16636 where timestamp between "2012-06-010 00:00" and "2013-07-02 23:59"  limit 10000000000000');

Temp_result = db.query();

Temp_result.timestamp(1);  % Temp & Pyr are almost syncronized - Thy have the same dates
Temp_result.value(1);

%from past to now
dt_Temp = datenum(Temp_result.timestamp(end:-1:1) );
Temp = Temp_result.value(end:-1:1) ;  %


% Brightness Left
%db.prepareStatement(['select * from id_16635 order by timestamp desc limit ',lim]); 
%db.prepareStatement('select * from id_16635 where timestamp between "2012-06-010 00:00" and "2013-07-02 23:59"  limit 1e8');

%BrightLeft_result = db.query();

% BrightLeft_result.timestamp(1); 
% BrightLeft_result.value(1) /1000;  % divide by 1000 to get Lux units

% Brightness Center
%db.prepareStatement(['select * from id_16634 order by timestamp desc limit ',lim]); 
%db.prepareStatement('select * from id_16634 where timestamp between "2012-06-010 00:00" and "2013-07-02 23:59"  limit 1e8');
%BrightCenter_result = db.query();

% BrightCenter_result.timestamp(1) 
% BrightCenter_result.value(1) /1000;  % divide by 1000 to get Lux units

% Brightness Right
%db.prepareStatement(['select * from id_16633 order by timestamp desc limit ',lim]); 
%db.prepareStatement('select * from id_16633 where timestamp between "2012-06-010 00:00" and "2013-07-02 23:59"  limit 1e8');
%BrightRight_result = db.query();

% BrightRight_result.timestamp(1); 
% BrightRight_result.value(1) /1000;  % divide by 1000 to get Lux units


% ----------------------  Power Generation Readings -----------------------
% Power
%db.prepareStatement(['select * from id_16525 order by timestamp desc limit ',lim]);
db.prepareStatement('select * from id_16525 where timestamp between "2012-06-010 00:00" and "2013-07-02 23:59"  limit 10000000000000');

PV_result = db.query();

%from past to now
dt_PV = datenum(PV_result.timestamp(end:-1:1) );
P_PV  = -1*PV_result.value(end:-1:1) ;  



%-------------- Plot Figures



figure(1)
%GHI
subplot(3,1,1)
plot(dt_GHI,GHI), 
%hold on, plot(dt_past_window,GHI_past,'m'), hold off 
title(['Measurement History of GHI for ', datestr(floor(dt_GHI(end)))])
ylabel('W/m^2')
xlim([dt_GHI(1) - 1/24/60, dt_GHI(end) + 1/24/60  ])
datetick('x','HH:MM','keeplimits')

% Temperature
subplot(3,1,2)
plot(dt_Temp,Temp),
ylabel ('�C')
title(['Measurement History of Temperature for ',datestr(floor(dt_GHI(end)))])
xlim([dt_GHI(1) - 1/24/60, dt_GHI(end) + 1/24/60  ])
datetick('x','HH:MM','keeplimits')

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




return







iter= 0;

%%  Get the last measurements


iter_stop = 1;
%while iter_stop == 1

iter= iter + 1

%------------- past windows




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
last_ten_min  = 10 /24/60:-1 /24/60:1 /24/60 ;

% Up to minuntes
now_time = day_now + hr_now/24 + min_now/24/60; 
%datestr(now_time) % Debug

% Past window
dt_past_window = now_time - last_ten_min';
%datestr(dt_past_window) % Debug








meas = 20; % Number of measurments required in window
lim = num2str (meas);

% ---------------------  Atmospheric Readings ----------------------------
%Weather Sensor Status
db.prepareStatement(['select * from id_16664 order by timestamp desc limit ',lim]);
WeatherSensorSatus_result2 = db.query();

%Weather  Status
db.prepareStatement(['select * from id_16660 order by timestamp desc limit ',lim]);
WeatherUnitStatus_result2 = db.query();

% Pyranometer
db.prepareStatement(['select * from id_16663 order by timestamp desc limit ',lim]);
%db.prepareStatement('select * from id_16663 where timestamp between "2013-06-10 00:00" and "2013-06-10 23:59"  limit 1000000');

GHI_result2 = db.query();


%from past to now
dt_GHI2 = datenum(GHI_result2.timestamp(end:-1:1) );
GHI2 = GHI_result2.value(end:-1:1)   * 100;  %Multiply by 100 to match KNX display 


% % %Wind Speed
db.prepareStatement(['select * from id_16637 order by timestamp desc limit ',lim]);
WindSpeed_result2 = db.query();


% % Wind Speed Threshold
% db.prepareStatement('select * from id_16641 order by timestamp desc limit ',lim]); 
% WindSpeedThreshold_result = db.query();


%Temperature
db.prepareStatement(['select * from id_16636 order by timestamp desc limit ',lim]); 
Temp_result2 = db.query();

%from past to now
dt_Temp2 = datenum(Temp_result2.timestamp(end:-1:1) );
Temp2 = Temp_result2.value(end:-1:1) ;  %


% % Brightness Left
% db.prepareStatement(['select * from id_16635 order by timestamp desc limit ',lim]); 
% BrightLeft_result2 = db.query();
% 
% % Brightness Center
% db.prepareStatement(['select * from id_16634 order by timestamp desc limit ',lim]); 
% BrightCenter_result2 = db.query();
% % divide values by 1000 to get Lux units
% 
% % Brightness Right
% db.prepareStatement(['select * from id_16633 order by timestamp desc limit ',lim]); 
% BrightRight_result2 = db.query();


% ----------------------  Power Generation Readings -----------------------
% Power
db.prepareStatement(['select * from id_16525 order by timestamp desc limit ',lim]);
PV_result2 = db.query();
 
%from past to now
dt_PV2 = datenum(PV_result2.timestamp(end:-1:1) );
P_PV2  = -1*PV_result2.value(end:-1:1)  ;  




%if iter > 1
    
%    if iter == 2
%         dt_GHI = [dt_GHI(1:end-meas); dt_GHI2];
%         GHI = [GHI(1:end-meas);   GHI2 ];   
% 
%         dt_Temp = [dt_Temp(1:end-meas); dt_Temp2 ];
%         Temp = [ Temp(1:end-meas);  Temp2 ];
% 
%         dt_PV = [dt_PV(1:end-meas); dt_PV2];
%         P_PV = [P_PV(1:end-meas); P_PV2];        
%    else
%         dt_GHI = [dt_GHI(1:end-meas); dt_GHI2];
%         GHI = [GHI(1:end-meas);   GHI2 ];   
% 
%         dt_Temp = [dt_Temp(1:end-meas); dt_Temp2 ];
%         Temp = [ Temp(1:end-meas);  Temp2 ];
% 
%         dt_PV = [dt_PV(1:end-meas); dt_PV2];
%         P_PV = [P_PV(1:end-meas); P_PV2];        
%         
%    end
    
    figure(1)
    %GHI
    subplot(3,1,1)
    hold on, plot(dt_GHI2,GHI2,'m'), hold off,
    %hold on, plot(dt_past_window,GHI_past,'m'), hold off 
    title(['Measurement History of GHI for ', datestr(floor(dt_GHI(end)))])
    ylabel('W/m^2')
    xlim([dt_GHI(1) - 1/24/60, dt_GHI(end) + 1/24/60  ])
    datetick('x','HH:MM','keeplimits')

    % Temperature
    subplot(3,1,2)
    hold on, plot(dt_Temp2,Temp2,'m'), hold off,
    ylabel ('�C')
    title(['Measurement History of Temperature for ',datestr(floor(dt_GHI(end)))])
    xlim([dt_GHI(1) - 1/24/60, dt_GHI(end) + 1/24/60  ])
    datetick('x','HH:MM','keeplimits')

    subplot(3,1,3)
    hold on, plot(dt_PV2,P_PV2,'m'), hold off,
    ylabel ('W')
    title(['Measurement History of Generated Power for ', datestr(floor(dt_PV(end)))])
    xlim([dt_GHI(1) - 1/24/60, dt_GHI(end) + 1/24/60])
    datetick('x','HH:MM','keeplimits')
    %set(gca,'XTick',0:0.2:1); %Set the XTicks values
    %datetick('x','HH','keepticks')
    %xlim([0 1])
    %ylim([0 1.2e3])
    
    
%end
















% Resample all variables using the past window
GHI_past  = interp1(dt_GHI2,GHI2,dt_past_window,'spline');
Temp_past = interp1(dt_Temp2,Temp2,dt_past_window,'spline');
P_PV_past = interp1(dt_PV2,P_PV2,dt_past_window,'spline');










figure(2)
subplot(3,2,1)
plot(dt_past_window,GHI_past)
xlim([dt_past_window(1) - 1/24/60, dt_past_window(end) + 1/24/60  ])
title(['10 min Past Window of GHI for ', datestr(floor(dt_GHI(end)))])
ylabel('W/m^2')
datetick('x','HH:MM','keeplimits')

% Temperature
subplot(3,2,3)
plot(dt_past_window,Temp_past),
xlim([dt_past_window(1) - 1/24/60, dt_past_window(end) + 1/24/60  ])
ylabel ('�C')
title(['10 min Past Window of Temperature for ',datestr(floor(dt_GHI(end)))])
datetick('x','HH:MM','keeplimits')

%Power
subplot(3,2,5)
plot(dt_past_window,P_PV_past), 
xlim([dt_past_window(1) - 1/24/60, dt_past_window(end) + 1/24/60  ])
ylabel ('W')
title(['10 min Past Window of Power for ', datestr(floor(dt_PV(end)))])
datetick('x','HH:MM','keeplimits')
%set(gca,'XTick',0:0.2:1); %Set the XTicks values
%datetick('x','HH','keepticks')
%xlim([0 1])
%ylim([0 1.2e3])




%end







% -------------- Apply Clear Sky Model --------------------------
Location.latitude =  8+16/60+38.50/60/60;     %   8�16'38.50" E   Local Longitude of CHCRC.C1 roof (negative if west)
Location.longitude = 47+27/60+33.92/60/60;    %   47�27'33.92" N  Local latitude of CHCRC.C1 roof  (positive if north of equator)
Location.altitude =  400; %455;                     %   [m]  


             
DN = dt_past_window; %GHIData{iter}.Time %+  start_date    ;

% Siwtzerland has a Central European Summer Time (CEST) +0200 UTC 
UTC = 2; %4.5; %  UTC is 2 but I have noticed an offset of 2:30 min.
Time = pvl_maketimestruct(DN, UTC);

% Clear sky GHI, DNI, and DHI from Ineichen/Perez model
[ClearSkyGHI_ineichen, ClearSkyDNI_ineichen, ClearSkyDHI_ineichen]= pvl_clearsky_ineichen(Time, Location);       

cloudyness_factor = 0.8 * mean(ClearSkyGHI_ineichen ./ GHI_past)  ;
modified_ClearSkyGHI_ineichen= ClearSkyGHI_ineichen ./ cloudyness_factor;

% Clear sky GHI from Haurwitz model
[SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(Time, Location);
ApparentZenith = 90-ApparentSunEl;

ClearSkyGHI_haurwitz = pvl_clearsky_haurwitz(ApparentZenith);
  
figure(2)
subplot(3,2,1)
hold on, plot(dt_past_window,modified_ClearSkyGHI_ineichen,'r'), hold off,
legend({'Measurement', 'Prediction'})
%xlim([dt_past_window(1) - 1/24/60, dt_past_window(end) + 1/24/60  ])









%------------- future windows

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
next_ten_min  = 1 /24/60:1 /24/60:10 /24/60 ;

% Up to minuntes
now_time = day_now + hr_now/24 + min_now/24/60; 
%datestr(now_time) % Debug

% Past window
dt_future_window = now_time + next_ten_min';
%datestr(dt_past_window) % Debug













% -------------- Apply Clear Sky Model for the Future --------------------------
Location.latitude =  8+16/60+38.50/60/60;     %   8�16'38.50" E   Local Longitude of CHCRC.C1 roof (negative if west)
Location.longitude = 47+27/60+33.92/60/60;    %   47�27'33.92" N  Local latitude of CHCRC.C1 roof  (positive if north of equator)
Location.altitude =  400; %455;                     %   [m]  


             
DN = dt_future_window; %GHIData{iter}.Time %+  start_date    ;

% Siwtzerland has a Central European Summer Time (CEST) +0200 UTC 
UTC = 2; %4.5; %  UTC is 2 but I have noticed an offset of 2:30 min.
Time = pvl_maketimestruct(DN, UTC);

% Clear sky GHI, DNI, and DHI from Ineichen/Perez model
[ClearSkyGHI_ineichen, ClearSkyDNI_ineichen, ClearSkyDHI_ineichen]= pvl_clearsky_ineichen(Time, Location);       

cloudyness_factor = 0.8 * mean(ClearSkyGHI_ineichen ./ GHI_past)  ;
modified_ClearSkyGHI_ineichen= ClearSkyGHI_ineichen ./ cloudyness_factor;

% Clear sky GHI from Haurwitz model
[SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(Time, Location);
ApparentZenith = 90-ApparentSunEl;

ClearSkyGHI_haurwitz = pvl_clearsky_haurwitz(ApparentZenith);








figure(2) % Use the same figure

subplot(3,2,2)
hold on, plot(dt_future_window,modified_ClearSkyGHI_ineichen,'r'),hold off
datetick('x','HH:MM','keeplimits')
xlim([dt_future_window(1) - 1/24/60, dt_future_window(end) + 1/24/60  ])
title(['10 min Future Window of GHI for ', datestr(floor(dt_GHI(end)))])
ylabel('W/m^2')


% Temperature
%subplot(3,2,4)
% plot(dt_past_window,Temp_past),
% xlim([dt_past_window(1) - 1/24/60, dt_past_window(end) + 1/24/60  ])
% ylabel ('�C')
% title(['10 min Past Window of Temperature for ',datestr(floor(dt_GHI(end)))])
% datetick('x','HH:MM','keeplimits')

%Power
%subplot(3,1,6)
% plot(dt_past_window,P_PV_past), 
% xlim([dt_past_window(1) - 1/24/60, dt_past_window(end) + 1/24/60  ])
% ylabel ('W')
% title(['10 min Past Window of Power for ', datestr(floor(dt_PV(end)))])
% datetick('x','HH:MM','keeplimits')
% %set(gca,'XTick',0:0.2:1); %Set the XTicks values
% %datetick('x','HH','keepticks')
% %xlim([0 1])
% %ylim([0 1.2e3])















%db.close();  %Remember to close the db when finish !!!!
