Location.latitude = 35.04;
Location.longitude = -106.62;
Location.altitude = 1619;
% Create 1-min time series for Jan 1, 2012
DN = datenum(2012, 6,1):1/(24*60):datenum(2012, 6, 1, 23, 59, 59);
Time = pvl_maketimestruct(DN, -7);
[ClearSkyGHI, ClearSkyDNI, ClearSkyDHI]= pvl_clearsky_ineichen(Time, Location);
dHr = Time.hour+Time.minute./60+Time.second./3600; % Calculate decimal hours for plotting
figure
plot(dHr,ClearSkyGHI)
hold all
plot(dHr,ClearSkyDNI)
plot(dHr,ClearSkyDHI)
title('Clear Sky Irradiance Components on June 1, 2012 in Albuquerque, NM')
xlabel('Hour of the Day (hr)')
ylabel('Irradiance (W/m^2)')
legend('GHI','DNI','DHI')