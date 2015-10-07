%% pvl_dirint 
% Determine DNI from GHI using the DIRINT modification of the DISC model
%
%% Syntax
%   |pvl_dirint = pvl_dirint(GHI, Z, doy, pressure)|
%   |pvl_dirint = pvl_dirint(GHI, Z, doy, pressure, UseDelKtPrime)|
%   |pvl_dirint = pvl_dirint(GHI, Z, doy, pressure, UseDelKtPrime,
%   DewPtTemp)|
%
%% Description
% Implements the modified DISC model known as "DIRINT" introduced in [1].
% DIRINT predicts direct normal irradiance (DNI) from measured global
% horizontal irradiance (GHI). DIRINT improves upon the DISC model by
% using time-series GHI data and dew point temperature information. The
% effectiveness of the DIRINT model improves with each piece of
% information provided.
%
%% Inputs:   
% * *|GHI|* - a scalar or vector of global horizontal irradiance in W/m^2. If GHI
%     is a vector it must be of the same size as all other vector inputs.
%     GHI must be >=0.
% * *|Z|* - a scalar or vector of true (not refraction-corrected) zenith
%     angles in decimal degrees. If Z is a vector it must be of the
%     same size as all other vector inputs. Z must be >=0 and <=180.
% * *|doy|* - a scalar or vector of values providing the day of the year. If
%     doy is a vector it must be of the same size as all other vector inputs.
%     doy must be >= 1 and < 367.
% * *|pressure|* - a scalar or vector of values providing the site pressure in
%     Pascal. If pressure is a vector it must be of the same size as all
%     other vector inputs. pressure must be >=0. Pressure may be measured
%     or an average pressure may be calculated from site altitude.
% * *|UseDelKtPrime|* - a numeric scalar indicating if the user would like to
%     utilize the time-series nature of the GHI measurements. A value of 0
%     will not use the time-series improvements, any other numeric value
%     will use time-series improvements. It is recommended that time-series
%     data only be used if the time between measured data points is less
%     than 1.5 hours. If |UseDelKtPrime| is not provided, the default is 1
%     (use time-series improvements). If none of the input arguments are
%     vectors, then time-series improvements are not used (because it's not
%     a time-series).
% * *|DewPtTemp|* - a scalar or vector of surface dew point temperatures, in 
%     degrees C. If DewPtTemp is a vector, it must be of the same size as
%     other vector inputs. Values of |DewPtTemp| may be numeric or NaN. Any
%     single time period point with a |DewPtTemp|=NaN does not have dew point
%     improvements applied. If |DewPtTemp| is not provided, then dew point 
%     improvements are not applied.  
%
%% Outputs:   
% * *|dirintDNI|* - the modeled direct normal irradiance in W/m^2 provided by the
%     DIRINT model. |dirintDNI| is a column vector with the same number of
%     elements as the input vector(s).
%
%% Example
% This example shows the measured and dirint-model-estimated direct normal
% irradiance for August 6 (TMY3 file for Albuquerque, NM).  Note that the
% model-estimated quantity is not a match.  This is a very difficult
% quantity to estimate and there are large uncertainties associated with
% using model estimates for DNI.
TMYData = pvl_readtmy3('723650TY.csv');
TimeMatlab = TMYData.DateNumber;
Time = pvl_maketimestruct(TimeMatlab, ones(size(TimeMatlab))*TMYData.SiteTimeZone); 
dayofyear = pvl_date2doy(Time.year, Time.month, Time.day);
DNI = TMYData.DNI; % Read in for comparison with results
DHI = TMYData.DHI; % Read in for comparison with results
GHI = TMYData.GHI;
Location = pvl_makelocationstruct(TMYData.SiteLatitude,TMYData.SiteLongitude,...
TMYData.SiteElevation); %Altitude is optional
PresPa = TMYData.Pressure*100; %Convert pressure from mbar to Pa
%
% Models that rely on sun position that are run at specific time steps
% (e.g., hourly) run into numerical problems during timesteps when the sun
% straddles the horizion (i.e., spends part of the timestep above and part
% below the horizion).  Since in this example we use TMY day, which is
% hourly and reported at the end of the hour, we adjust the sun position in
% the following two ways.  (1) For hours when the sun is above the horizion
% we adjust sun position so it is for the midle of the hour.  For example,
% sun position reported at 4PM is the sun position at 3:30 PM.  (2) For
% hours where the sun straddles the horizion, we report the position half
% way between the horizion and the position of the sun at the end of the
% hour.  This is appropriate because the GHI data is essentially a sum of the GHI
% measurements made in the past hour.
%
%Run sun position twice:  Once for end of hour positions.  
% A second time for mid hour positions.  
% Adjust sun elevation and GHI for hours when the sun traverses the horizion
% Include these adjusted values in the mid hour positions
%
[~, ~, AppSunEl, ~] = pvl_ephemeris(Time,Location,PresPa,TMYData.DryBulb);
Time.hour = Time.hour-.5; % shift times back 1/2 hour for sun position calculation because of tmy
% timestamps
[~, ~, AppSunEla, ~] = pvl_ephemeris(Time,Location,PresPa,TMYData.DryBulb);

A=diff(sign(AppSunEl)); %identifies hour before sun straddles horizion (2,-2)
ind1 = find(A==2)+1; % AM hour where sun straddles horizion
ind2 = find(A==-2)+1; % PM hour where sun straddles horizion
%AM Adjustment
    AppSunEl(ind1) = AppSunEl(ind1)/2; %change sun elevation to mid way above the horizion 
 %PM Adjustment
    AppSunEl(ind2) = AppSunEl(ind2)/2; %change sun elevation to mid way above the horizion 
AppSunEla(ind1)=  AppSunEl(ind1);
AppSunEla(ind2)=  AppSunEl(ind2);
DNI_model = pvl_dirint(GHI,90-AppSunEla, dayofyear, PresPa);

figure
tfilter = and(Time.month == 8,Time.day == 6);
plot(Time.hour(tfilter),DNI_model(tfilter),'-s')
hold all
plot(Time.hour(tfilter),DNI(tfilter),'-o')
legend('DNI (DIRINT Model)','DNI (TMY3)','Location','NW')
xlabel('Hour of Day')
ylabel('Irradiance (W/m^2)')
title('Albuquerque Direct Normal Irradiance Comparison - Aug 6','FontSize',14)
%% References
%
% [1]Perez, R., P. Ineichen, E. Maxwell, R. Seals and A. Zelenka, (1992).
%   "Dynamic Global-to-Direct Irradiance Conversion Models".  ASHRAE 
%   Transactions-Research Series, pp. 354-369
%
% [2] Maxwell, E. L., "A Quasi-Physical Model for Converting Hourly 
%   Global Horizontal to Direct Normal Insolation", Technical 
%   Report No. SERI/TR-215-3087, Golden, CO: Solar Energy Research 
%   Institute, 1987.
%
% DIRINT model requires time series data (ie. one of the inputs must be a
% vector of length >2.
%
%% See also
% <pvl_disc_help.html |pvl_disc|> , 
% <pvl_erbs_help.html |pvl_erbs|> , 
% <pvl_louche_help.html |pvl_louche|> , 
% <pvl_orgill_hollands_help.html |pvl_orgill_hollands|> , 
% <pvl_reindl_1_help.html |pvl_reindl_1|> , 
% <pvl_reindl_2_help.html |pvl_reindl_2|> , 
% <pvl_ephemeris_help.html |pvl_ephemeris|> ,
% <pvl_date2doy_help.html |pvl_date2doy|> ,
% <pvl_alt2pres_help.html |pvl_alt2pres|> 

%%
% Copyright 2014 Sandia National Laboratories

