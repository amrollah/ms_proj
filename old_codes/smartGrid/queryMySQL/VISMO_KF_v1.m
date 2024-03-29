%% Main File: VISMO_KF_v1.m
%
% A. Cortinovis, andrea.cortinovis@ch.abb.com
% ABB CHCRC, 12.08.2013, Project: VISMO
%
% Summary: EKF implementation using past data from the smartgrid SQL-database 
% and predction for the next 15min
%
% Dependencies: one_KF_step.m, queryMYSQL library, WLAN connection to smartgrid 
%               demo lab
%
% -------------------------------------------------------------------------
% clear all
% clear functions
close all
% clf
clc
test = false;
queryDatabase = true;
getPVandGHIFromMat = false;
% getPredictedIrrValuesFromMat = true;
if(getPVandGHIFromMat)
    addpath('../../VISMO_Code'); %Adapt the directory for your own directories.
    addpath(genpath('../../VISMO_Code/Utils/Data'));
    load 'PVGHIData_plus_CSK.mat';
    dayOfTheMonth = 12; %we're interested with August 12th 2013.
    %get the Irradiance value prediction from the VISMO
    %multiObjectTrackingCV algorithm.
%     load 'estimatedIrradianceValues.mat'; %don't forget to change it with ~temporary version.
    load 'estimatedIrradianceValuesTemporary2.mat'; 
    shiftTime = 0;
    if(test) % test 5 mins of shift in timestamp.
        shiftTime = 5/60/24;
    end
    for i=1:length(estimatedIrradianceValues)
        estimatedIrradianceValues(i).timeStamp = datenum(estimatedIrradianceValues(i).dateStr) - shiftTime;
    end
    g = GHIData{dayOfTheMonth};
    GHI_result.timestamp = datestr(g.Time + g.DateNumber);
    dt_GHI = g.Time + g.DateNumber;
%     dt_GHI = dt_GHI; %reverse the order as done with the smart grid database
    GHI_result.value = g.GHI; %notice that this is already multiplied by 100 to match the KNX display
    GHI = GHI_result.value;
    PV_result.timestamp = GHI_result.timestamp; %note they are all interpolated in the saved .mat file.
    if(1) 
        PV_result.value = g.PV_output/1000; %scaling the values to kW from W. 
    else %we are scaling the power.
        PV_result.value = g.PV_output;
    end
    dt_PV = dt_GHI;
    P_PV = PV_result.value;
    
    dt_Temp = dt_GHI;
    Temp = g.ambTemp;
    totalTimeLapse = datevec(  abs(dt_GHI(end)-dt_GHI(1))); %taking only up until several days of length into consideration. for longer time queries, please modify the code.
    get_hours = totalTimeLapse(3)*24 + totalTimeLapse(4) + totalTimeLapse(5)/60 + totalTimeLapse(6)/60/60;
end
% load libs
addpath(genpath('../PV_LIB Version 1_1'));
%% Set up SQL connection
addpath(fullfile(pwd, 'src'));


%% Get Data from SQL db
if(queryDatabase)
    javaaddpath('lib/mysql-connector-java-5.1.6/mysql-connector-java-5.1.6-bin.jar');
        % import classes
    import edu.stanford.covert.db.MySQLDatabase;

        % create database connection
    db = MySQLDatabase('10.41.94.39', 'smartgrid', 'root', 'smartgrid');
    % don't forget to close the db at the end!

    tic;
    get_hours = 1/6;% amount of hours you want to get from database.

    % Pyranometer (internal sampling rate is 5sec)
    lim_GHI = get_hours; % in hours
    lim = num2str (round(lim_GHI * 60 * 60 / 5) );

%     dbTimeInterval = ' where timestamp between "2014-08-01 05:01" and "2014-09-01 15:30"  limit 1000000';

    db.prepareStatement(['select * from id_16663 order by timestamp desc limit ',lim]);
%     db.prepareStatement(['select * from id_16663' dbTimeInterval]);
    GHI_result = db.query();
    dt_GHI = datenum(GHI_result.timestamp(end:-1:1) ); %from past to now
    GHI = GHI_result.value(end:-1:1)  * 100;  %Multiply by 100 to match KNX display 


    % % % % % %Wind Speed (internal sampling rate is 5-6sec) 
    % % % db.prepareStatement(['select * from id_16637 order by timestamp desc limit ',lim]);
    % % % WindSpeed_result = db.query(); % Wind does not work !!
    % % % Wind Speed Threshold
    % % db.prepareStatement('select * from id_16641 order by timestamp desc limit ',lim]); 
    % % WindSpeedThreshold_result = db.query();
    % %  
    % % WindSpeedThreshold_result.timestamp(1)  
    % % WindSpeedThreshold_result.value(1)

    % Temperature (internal sampling rate is 5sec)
    lim_Temp = get_hours; % in hours
    lim = num2str (round(lim_Temp * 60 * 60 / 5) );
%     db.prepareStatement(['select * from id_16636' dbTimeInterval]); 
    db.prepareStatement(['select * from id_16636 order by timestamp desc limit ',lim]);
    Temp_result = db.query();
    dt_Temp = datenum(Temp_result.timestamp(end:-1:1) );%from past to now
    Temp = Temp_result.value(end:-1:1) ;  %from past to now


    % Power (internal sampling rate is 80sec)
    lim_PV = get_hours; % in hours
    lim = num2str (round(lim_PV * 60 * 60 / 80) );
%     db.prepareStatement(['select * from id_16525' dbTimeInterval]);
    db.prepareStatement(['select * from id_16525 order by timestamp desc limit ',lim]);
    PV_result = db.query();
    dt_PV = datenum(PV_result.timestamp(end:-1:1) ); %from past to now
    P_PV = -1 * PV_result.value(end:-1:1) ;  

    time = toc;
    disp(['Info: Data Acquisition time: ' num2str(time) ' sec'])
    
    db.close();  %Remember to close the db when finish !!!!
end
%%create the Accuracy variables
% accuracyDynamicWindow.MBE = {};
% accuracyDynamicWindow.RMSE = {};
accuracyFixedWindow = {};
accuracyFixedWindow.DistanceIntoFuture = {};
accuracyFixedWindow.RMSE = {};
accuracyFixedWindow.MBE = {};
accuracyFixedWindow.MAE = {};

%% Create Time series
tic; 

GHI_ts = timeseries(GHI,datestr(dt_GHI));
set(GHI_ts,'Name','GHI in [W/m2]');

PV_ts = timeseries(P_PV,datestr(dt_PV));
set(PV_ts,'Name','Power in [W]');

Temp_ts = timeseries(Temp,datestr(dt_Temp));
set(Temp_ts,'Name','Temp in [�C]');

%% Plot Figures (unsynchronized, in blue)

% figure(1)
% subplot(3,1,1)
% plot(GHI_ts) 
% subplot(3,1,2)
% plot(Temp_ts)
% subplot(3,1,3)
% plot(PV_ts)

%% Loop over the data to have chunks of 10 mins
% Check the earliest possible time we have from the Algorithm. (hopefully
% available from the start).
timeStampIrradiance = zeros(length(estimatedIrradianceValues),1);
for i=1:length(estimatedIrradianceValues)
   timeStampIrradiance(i) = estimatedIrradianceValues(i).timeStamp;
end
%note that here I'm assuming that GHI from sensors have wider time interval
%of recording (i.e. GHI_ts is from 24/7 data but timeStampIrradiance is
%7/7)

%making it start around 8.10..
[dist, startingIndex] = min(abs(timeStampIrradiance(50)-dt_GHI));
display(['Method time: ' datestr(timeStampIrradiance(50)) ', GHI time from smartgrid: ' datestr(dt_GHI(startingIndex)) ]);
%notice I'm changing get_hours again to fit to 10 mins lapse.
get_hours = 1/6; %10 mins
%check if the sensors have the starting time -10 mins available.
timeNeededFromPast = datenum([0, 0, 0, floor(get_hours), floor(mod(get_hours,1)*60), mod(mod(get_hours,1)*60,1)*60]);
if(min(dt_GHI) > dt_GHI(startingIndex) - timeNeededFromPast) %we need to shift time since sensors don't have sufficient 'PAST' information
    startingIndex = length(dt_GHI); %note that dt_GHI is reversed.
end
stillPast = true;
it = 1;
clear power5MinsLater IrradiationValues IrrEstimated;
power5MinsLater = {};
next_15_min  = 1/24/60/60 : 1/24/60/60*5 : 15/24/60;
next_5_min  = 1/24/60/60 : 1/24/60/60*5 : 5/24/60; %every 5 seconds


%%here, 'ind' points to starting index - 10 mins. So 10 mins into the past.
% [val, ind] = min(abs(dt_GHI - (dt_GHI(startingIndex) -
% timeNeededFromPast))); % start the recording from earliest possible.
ind = 1;
while stillPast
    power5MinsLater{it}.Time = datestr(dt_GHI(ind+it-1));
    power5MinsLater{it}.MeasuredPower = P_PV(ind+it-1); %note that power ghi and all have the same timestamps as they are interpolated.
    power5MinsLater{it}.EstimatedPower = nan; %TODO: Notice for now I'm using 0, later nan!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IrradiationValues(it) = GHI_ts.Data(ind+it-1);
    IrrEstimated(it) = nan; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(ind+it == startingIndex)
       stillPast = false; 
    end
    it = it+1;
end

%we currently only keep current measurements for future. Though this can be
%improved. Only for the initial 5 mins though, so doesn't matter that much.
stillLessThan5MinIntoFuture  = true;
ind = startingIndex;
%Here, 'ind' points to starting index, so this very moment. Loop continues
%until Now + (5mins - 1frame) into future.
while stillLessThan5MinIntoFuture
    power5MinsLater{end+1}.Time = datestr(dt_GHI(ind));
    power5MinsLater{end}.MeasuredPower = P_PV(ind);
    power5MinsLater{end}.EstimatedPower = nan;%TODO: Notice for now I'm using the measured power instead! 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    IrradiationValues(end+1) = GHI_ts.Data(ind);
    IrrEstimated(end+1) = GHI_ts.Data(ind);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(dt_GHI(startingIndex) + next_5_min(end) <= dt_GHI(ind-1))
        stillLessThan5MinIntoFuture = false;
        break;
    end
    ind = ind+1;
end

finishLoop = false;
iterationIndex = 1;
loopStartingIndex = startingIndex;
while ~finishLoop %return here FIX ME.
    now_vec = datevec(dt_GHI(startingIndex));
    
    %% Resampling and synchonization
    
    % now_vec = datevec(dt_GHI(1)); %was previously now_vec = datevec(dt_GHI(end));
    
    day_now = datenum(now_vec(1:3));
    %datestr(day_now) % debug
    hr_now  =  (now_vec(4));
    min_now = (now_vec(5)) ;
    sec_now = (now_vec(6));
    
    %Don't get confused with var name, it's basically time into the past.
    % "get_hours-1" interval 1min resolution in units of 'days' with
    % sampling rate of 5 seconds
    last_ten_min  = (get_hours)*60/24/60 : -1/24/60/60*5 : 0/24/60; %ASK WHY it was get_hours-1!
    
    % Up to minuntes
    now_time = day_now + hr_now/24 + min_now/24/60 + sec_now/24/60/60;
    %datestr(now_time) % Debug
    
    % Past window
    dt_past_window = now_time - last_ten_min';
    %datestr(dt_past_window) % Debug
    
    % Resample
    GHI_past  = interp1(dt_GHI,GHI,dt_past_window,'spline');
    Temp_past = interp1(dt_Temp,Temp,dt_past_window,'spline');
    P_PV_past = interp1(dt_PV,P_PV,dt_past_window,'spline');
    
    % Create synchronized Time series
    GHI_ts2 = timeseries(GHI_past,datestr(dt_past_window));
    set(GHI_ts2,'Name','GHI in [W/m2]');
    
    PV_ts2 = timeseries(P_PV_past,datestr(dt_past_window));
    set(PV_ts2,'Name','Power in [W]');
    
    Temp_ts2 = timeseries(Temp_past,datestr(dt_past_window));
    set(Temp_ts2,'Name','Temp in [�C]');
    
    %% Plot in same figure (synchronized, in magenta)
%     figure(1)
%     subplot(3,1,1)
%     hold on, plot(GHI_ts2,'color','m')
%     
%     subplot(3,1,2)
%     hold on, plot(Temp_ts2,'m')
%     ylabel('Temp in [�C]')
%     
%     subplot(3,1,3)
%     hold on, plot(PV_ts2,'m')
    
    %% Generate Randomized prediction for next 15min
    end_point = GHI_ts2.Data(end);
    amplitude = end_point/10;
    
    next_15_min  = 1/24/60/60*5 : 1/24/60/60*5 : 15/24/60;
    next_5_min  = 1/24/60/60*5 : 1/24/60/60*5 : 5/24/60; %every 5 seconds
    now_time = day_now + hr_now/24 + min_now/24/60 + sec_now/24/60/60;
    dt_next_window = now_time + next_5_min';
    
    
    %PV measurements of future 
    if(dt_GHI(startingIndex) + next_5_min(end) <= dt_GHI(end))
        PVFuture= interp1(dt_PV,P_PV,dt_next_window,'spline');
        PV_futureMeasurements = timeseries(PVFuture,datestr(dt_next_window));
    end
    
    scenario = 'constant';
    switch scenario
        case 'falling'
            IRR_pred = end_point + [- [[1:15]*amplitude]' + randn(15,1)];
        case 'raising'
            IRR_pred = end_point + [ [[1:15]*amplitude]' + randn(15,1)];
        case 'constant'
            IRR_pred = end_point + randn(length(dt_next_window),1);
        case 'custom'
            %         IRR_pred = [end_point 890 1000 600 300 400 300 200 300 1000 1100 990 993 1000 1100];
            IRR_pred = [end_point 400 401 350 300 270 265 270 272 271 300 350 401 403 402];
        case 'clearSky'
            % add paths:
            str1 = 'X:\Projects\VISMO\PV_LIB Version 1_1';
            str2 = 'X:\Projects\VISMO\PV_LIB Version 1_1\Required Data';
            addpath(str1)
            addpath(str2)
            
            % -------------- Apply Clear Sky Model --------------------------
            Location.latitude =  47+27/60+33.92/60/60;    %   47�27'33.92" N  Local latitude of CHCRC.C1 roof  (positive if north of equator)
            Location.longitude = 8+16/60+38.50/60/60;     %   8�16'38.50" E   Local Longitude of CHCRC.C1 roof (negative if west)
            Location.altitude =  455;                     %   [m]
            
            DN = dt_next_window; %GHIData{iter}.Time %+  start_date    ;
            UTC = 2; %4.5; %  UTC is 2 but I have noticed an offset of 2:30 min.
            Time = pvl_maketimestruct(DN, UTC);
            
            % Clear sky GHI, DNI, and DHI from Ineichen/Perez model
            [ClearSkyGHI_ineichen, ClearSkyDNI_ineichen, ClearSkyDHI_ineichen]= pvl_clearsky_ineichen(Time, Location);
            N_temp = length(ClearSkyGHI_ineichen);
            temp_mean = ClearSkyGHI_ineichen ./ GHI_ts2.Data(end-N_temp+1:end);
            index2 = ~isnan(temp_mean);
            cloudyness_factor = 0.8 * mean(temp_mean(index2))  ;
            modified_ClearSkyGHI_ineichen= ClearSkyGHI_ineichen ./ cloudyness_factor;
            IRR_pred = modified_ClearSkyGHI_ineichen*367/428.8;
        otherwise
            error('No valid scenario is selected');
    end
    
    
    
    pred_ts = timeseries(IRR_pred,datestr(dt_next_window));
%     set(pred_ts,'Name','GHI in [W/m2]');
%     figure(1)
%     subplot(3,1,1)
%     hold on;
%     plot(pred_ts,':r')
%     ylabel('GHI in [W/m2]')
    
    time = toc;
%     disp(['Info: Object Creation and Plotting time: ' num2str(time) ' sec'])
    
    %% KF implementation last N samples and over next M samples
    
    tic;
    N = length(dt_past_window);
    M = length(dt_next_window);
    TEMP = zeros(N+M,1);
    TEMP2 = zeros(N+M,2);
    TEMP3 = zeros(N+M,2);
    [~, inde] = min(abs(timeStampIrradiance - now_time)); 
    %loop over past 10 mins until 5mins future.
    for i=1:N+M % loop over all samples
        if i <=N
            % historical data
            %train KF with past 10 mins.
            z_data = [GHI_ts2.Data(i) Temp_ts2.Data(i) PV_ts2.Data(i)];
        else
            % predicted data
%             z_data = [pred_ts.Data(i-N) NaN NaN];
            %predict future 5 mins using the trained KF.
%          %%%%%%%%%%%% %%%%%%%%%%%%%
%             [~, inde] = min(abs(timeStampIrradiance - dt_next_window(i-N))); 
            pos = (i-N)/M;
            totLenPrediction = length(estimatedIrradianceValues(inde).values);
            futureElement = ceil(min(pos,1)*totLenPrediction);
            GHI_est = estimatedIrradianceValues(inde).values(futureElement); 
            z_data = [GHI_est nan nan];
            timeStampCurrent = estimatedIrradianceValues(inde).timeStamp;
%          %%%%%%%%%%%%%
        end
        
        [x_new, y, K] = one_KF_step(z_data);
        if(y>5)
            y=5;
        elseif(y<0)
            y=0;
        end
        TEMP(i,1) = y;
        TEMP2(i,1:2) = x_new;
        TEMP3(i,1:2) = K;
        %NOTE that in order for the if below to work, you need to have
        %dt_next_window <= 5 mins
%         if(i > N-1 && getPredictedIrrValuesFromMat) %fill in the rest of Power values with the Power formula and predicted Irradiance values.
%             for k=N+1:N+M
%                 TEMP3(k,1:2) = K; %Not changing.
%                 TEMP2(k,1:2) = x_new; %Not changing.
%                 %find the corresponding GHI from future among
%                 %estimatedIrradianceValues
%                 
%                 %USE dt_next_window(k-N) for the timestamp of corresponding future timepoint,
%                 %compare with estimatedIrr....timeStamp closest match with
%                 %min fct's index. then maybe even interpolate the estimation if
%                 %you want high precision.
%                 
%                 [distance, index] = min(abs(timeStampIrradiance - dt_next_window(k-N)));
% 
%                 GHI_t = estimatedIrradianceValues(index).values(end); % Note that currently estimatedIrradianceValues only cover information for the next 5 mins where last element is the 5th minute prediction
%                 timeStampCurrent = estimatedIrradianceValues(index).timeStamp;
%                 TEMP(k,1) = x_new(1)*GHI_t + x_new(2)*(GHI_t)^2;
%             end
%             
%             break;
%         end
    end
    
    %Edit the variable that will contain the Daily History
    power5MinsLater{end+1}.Time = datestr(timeStampCurrent); %format is 12-Aug-2013 11:29:20 to prevent confusion. Time of 5 mins later
    power5MinsLater{end}.EstimatedPower = TEMP(end); %estimated power at the end of estimation (currently 5 mins) Estimated power corresponding to 5 mins later.
    
    [~, ind] = min(abs(timeStampCurrent-dt_PV)); % should correspond to dt_next_window(end)
    power5MinsLater{end}.MeasuredPower = P_PV(ind); %measured power of 5 mins later.
    IrradiationValues(end+1) = GHI_ts.Data(ind); %measured irradiance
    IrrEstimated(end+1) = GHI_est; %estimated irradiance value.
    
    
    
    
    %%%%%%%% SLOWS DOWN THE SPEED, DON'T DO IT FOR EACH DEBUG %%%%%%
    %Create the tsv file that is for the next 5 mins only.
    [~, ind] = min(abs(dt_past_window(1) - timeStampIrradiance));
%     fileName = ['../powerPrecise/Power_' datestr(timeStampIrradiance(ind),'yyyy-mm-dd_HH-MM-SS') '.tsv']; %Time is -10 mins to the past.
    fileName = ['C:/Users/CHFIOEZ/programs/xampp/htdocs/plotData/powerPrecise/Power_' datestr(dt_past_window(1),'yyyy-mm-dd_HH-MM-SS') '.tsv']; %Time is -10 mins to the past.
    clear PowerPastWindowPlus5MinsFuture;
    PowerPastWindowPlus5MinsFuture.Time = [dt_past_window; dt_next_window];
    PowerPastWindowPlus5MinsFuture.MeasuredPower = PV_ts2.Data;
    PowerPastWindowPlus5MinsFuture.PredictedPower = TEMP;
%     
    colNames = {{'times'}, {'Measured Power'}, {'Predicted Power'}};
    saveDatatoTSVPrecise(PowerPastWindowPlus5MinsFuture, colNames, fileName);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%Analysis on the RMSE and MBE for future estimations of PV Power%%%%
    if(dt_GHI(startingIndex)  && dt_GHI(startingIndex) + next_5_min(end) <= dt_GHI(end))
        accuracyDynamicWindow.RMSE(iterationIndex) = RMSE(PowerPastWindowPlus5MinsFuture.PredictedPower(N+1:end), PV_futureMeasurements.Data); %notice N = length(dt_past_window)
        accuracyDynamicWindow.MBE(iterationIndex) = MBE(PowerPastWindowPlus5MinsFuture.PredictedPower(N+1:end), PV_futureMeasurements.Data);
        accuracyDynamicWindow.MAE(iterationIndex) = MAE(PowerPastWindowPlus5MinsFuture.PredictedPower(N+1:end), PV_futureMeasurements.Data);
        %x Min ahead Prediction accuracy
        maxFuturePoint = 5; %this is a constant set for next 5 min.
        windowSize = ceil(1/maxFuturePoint*length(dt_next_window));
        halfWS = ceil(windowSize/3);
%         ind = length(accuracyFixedWindow.DistanceIntoFuture{1});
        for x = 1:5
            interestedFuturePoint = x;
            index = ceil(interestedFuturePoint/maxFuturePoint*length(dt_next_window));
            accuracyFixedWindow.DistanceIntoFuture{x}= x; %TODO: this is very non-necessary repetition. put this out of the loop.
            if(x ~= 5)
                accuracyFixedWindow.RMSE{x}(iterationIndex) = RMSE(PowerPastWindowPlus5MinsFuture.PredictedPower(N+index-halfWS:N+index+halfWS), PV_futureMeasurements.Data(index-halfWS:index+halfWS));
                accuracyFixedWindow.MBE{x}(iterationIndex) = MBE(PowerPastWindowPlus5MinsFuture.PredictedPower(N+index-halfWS:N+index+halfWS), PV_futureMeasurements.Data(index-halfWS:index+halfWS));
                accuracyFixedWindow.MAE{x}(iterationIndex) = MAE(PowerPastWindowPlus5MinsFuture.PredictedPower(N+index-halfWS:N+index+halfWS), PV_futureMeasurements.Data(index-halfWS:index+halfWS));
            else %we can't center 5th minute
                accuracyFixedWindow.RMSE{x}(iterationIndex) = RMSE(PowerPastWindowPlus5MinsFuture.PredictedPower(N+index-halfWS:N+index), PV_futureMeasurements.Data(index-halfWS:index));
                accuracyFixedWindow.MBE{x}(iterationIndex) = MBE(PowerPastWindowPlus5MinsFuture.PredictedPower(N+index-halfWS:N+index), PV_futureMeasurements.Data(index-halfWS:index));
                accuracyFixedWindow.MAE{x}(iterationIndex) = MAE(PowerPastWindowPlus5MinsFuture.PredictedPower(N+index-halfWS:N+index), PV_futureMeasurements.Data(index-halfWS:index));
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % %% Plot KF Restuls - Corrected + Future Prediction
    % figure;
    % subplot(3,1,1)
    % plot(TEMP)
    % hold on;
    % plot(N+1:M+N,TEMP(N+1:M+N),'r')
    % ylabel('PV Output [W]')
    % subplot(3,1,2)
    % plot(TEMP2(:,1))
    % hold on;
    % plot(N+1:M+N,TEMP2(N+1:M+N,1),'r')
    % ylabel('KF state \alpha_1')
    % subplot(3,1,3)
    % plot(TEMP2(:,2))
    % hold on;
    % plot(N+1:M+N,TEMP2(N+1:M+N,2),'r')
    % ylabel('KF state \alpha_2')
    %
    % KF_past_ts = timeseries(TEMP(1:N),datestr(dt_past_window));
    % set(KF_past_ts,'Name','PV_out in [W]');
    % figure(1)
    % subplot(3,1,3)
    % hold on;
    % plot(KF_past_ts,'r')
    %
    % KF_next_ts = timeseries(TEMP(N+1:M+N),datestr(dt_next_window));
    % set(KF_next_ts,'Name','PV_out in [W]');
    % figure(1)
    % subplot(3,1,3)
    % hold on;
    % plot(KF_next_ts,':r')
    % legend('Raw Signal','Synchonized Signal','KF Estimation','KF 15min Forecast','Location','northwest')
    % ylabel('Power in [W]')
    
    time = toc;
%     disp(['Info: KF execution and Plotting time: ' num2str(time) ' sec'])
    
    if(timeStampCurrent == estimatedIrradianceValues(end).timeStamp) %if we came to the end of predicted Irradiance values into the future, exit loop. Power estimation into future is over.
        finishLoop = true;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     THIS IF IS TEMPORARY. FOR DEBUGGING ONLY. REMOVE.
%     if(startingIndex == 14000) %if we came to the end of predicted Irradiance values into the future, exit loop. Power estimation into future is over.
%         finishLoop = true;
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    startingIndex = startingIndex+1;
    iterationIndex  = iterationIndex +1;
    display(['Starting Index of GHI_ts2 is ' num2str(startingIndex)]);
end
%save only the 5th min. estimation of each sliding window.
% fileName = 'C:/Users/CHFIOEZ/programs/xampp/htdocs/plotData/power.tsv';
fileName = '../powerTEST.tsv';
colNames = {{'times'}, {'Measured Power'}, {'Predicted Power'}};
% saveDatatoTSV(power5MinsLater, colNames, fileName);

clear mPower ePower timePower;
for i = 1:length(power5MinsLater)
    mPower(i) = power5MinsLater{i}.MeasuredPower;
    ePower(i) = power5MinsLater{i}.EstimatedPower;
    if((mPower(i)<0))
        ePower(i) = -.003;
    elseif(isnan(ePower(i)))
        ePower(i) = 0;
    end
    timePower(i) = datenum(power5MinsLater{i}.Time);
end
mPow_ts = timeseries(mPower,datestr(timePower));
ePow_ts = timeseries(ePower,datestr(timePower));
irr_ts = timeseries(IrradiationValues, datestr(timePower));
irrEst_ts = timeseries(IrrEstimated, datestr(timePower));
figure(1), subplot(2,1,1),
plot(mPow_ts, '.:r'); legend('Measured Power'); ylim([0, 5]);
subplot(2,1,2),
plot(ePow_ts, '.:b'); legend('Estimated Power'); ylim([0, 5]);
figure(2)
subplot(2,1,1),
plot(irr_ts, ':.m');  legend('Irradiation Values')
subplot(2,1,2),
plot(irrEst_ts, '.:r');  legend('Estimated Irradiation Values')

%%Plot RMSE and MBE
dynamicRMSE = timeseries(accuracyDynamicWindow.RMSE, datestr(dt_GHI(loopStartingIndex:startingIndex-1)));
figure(3)
plot(dynamicRMSE, ':.'); title('RMSE of dynamic windows');
figure(4), title('RMSE of fixed windows');
for i = 1:5
    subplot(5,1,i),
    fixedRMSE{i} = timeseries(accuracyFixedWindow.RMSE{i}, datestr(dt_GHI(loopStartingIndex:startingIndex-1)));
    plot(fixedRMSE{i}, ':.');
   %Note that I'm not plotting MBE yet. 
end

dynamicMAE = timeseries(accuracyDynamicWindow.MAE, datestr(dt_GHI(loopStartingIndex:startingIndex-1)));
figure(5)
plot(dynamicMAE, ':.'); title('MAE of dynamic windows');

%%find confidence interval:
[sortedDynamicRMSE indices]= sort(accuracyDynamicWindow.RMSE);
interestedInterval = ceil(length(sortedDynamicRMSE)*.95);
confDynamicRMSE = mean(sortedDynamicRMSE(1:interestedInterval))

display(sprintf (['Accuracy tests complete. Results are as follows:\nDynamic Window: \n NRMSE: ' num2str(mean(accuracyDynamicWindow.RMSE)) ',\t MBE: ' num2str(mean(accuracyDynamicWindow.MBE)) '\n Fixed Window 1 min forecast:\n RMSE: ' ...
    num2str(mean(accuracyFixedWindow.RMSE{1})) ',\t MBE: ' num2str(mean(accuracyFixedWindow.MBE{1})) '\n Fixed Window 2 min forecast:\n RMSE: ' num2str(mean(accuracyFixedWindow.RMSE{2})) ',\t MBE: ' ...
    num2str(mean(accuracyFixedWindow.MBE{2})) '\n Fixed Window 3 min forecast:\n RMSE: ' num2str(mean(accuracyFixedWindow.RMSE{3})) ',\t MBE: ' num2str(mean(accuracyFixedWindow.MBE{3})) ...
    '\n Fixed Window 4 min forecast:\n RMSE: ' num2str(mean(accuracyFixedWindow.RMSE{4})) ',\t MBE: ' num2str(mean(accuracyFixedWindow.MBE{4})) ...
    '\n Fixed Window 5 min forecast:\n RMSE: ' num2str(mean(accuracyFixedWindow.RMSE{5})) ',\t MBE: ' num2str(mean(accuracyFixedWindow.MBE{5})) ]));

% End of Code

%%
% 
return





%%

% str1 =  'X:\Projects\VISMO\PV_LIB Version 1_1';
% addpath(str1);
% -------------- Apply Clear Sky Model --------------------------
Location.latitude =  47+27/60+33.92/60/60;    %   47�27'33.92" N  Local latitude of CHCRC.C1 roof  (positive if north of equator)
Location.longitude = 8+16/60+38.50/60/60;     %   8�16'38.50" E   Local Longitude of CHCRC.C1 roof (negative if west)
Location.altitude =  455;                     %   [m]  
             
DN = dt_next_window; %GHIData{iter}.Time %+  start_date    ;

% Siwtzerland has a Central European Summer Time (CEST) +0200 UTC 
UTC = 2; %4.5; %  UTC is 2 but I have noticed an offset of 2:30 min.
Time = pvl_maketimestruct(DN, UTC);

% Clear sky GHI, DNI, and DHI from Ineichen/Perez model
[ClearSkyGHI_ineichen, ClearSkyDNI_ineichen, ClearSkyDHI_ineichen]= pvl_clearsky_ineichen(Time, Location);       

% index = find(GHI_ts2.Data<1);
% pred_ts.Data(index) = NaN;
N_temp = length(ClearSkyGHI_ineichen);
temp_mean = ClearSkyGHI_ineichen ./ GHI_ts2.Data(end-N_temp+1:end);
index2 = ~isnan(temp_mean);
cloudyness_factor = 0.8 * mean(temp_mean(index2))  ;
modified_ClearSkyGHI_ineichen= ClearSkyGHI_ineichen ./ cloudyness_factor;
% 
% % Clear sky GHI from Haurwitz model
[SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(Time, Location);
ApparentZenith = 90-ApparentSunEl;

ClearSkyGHI_haurwitz = pvl_clearsky_haurwitz(ApparentZenith);
  


% CS_pred_ts = timeseries(dt_next_window,modified_ClearSkyGHI_ineichen);
CS_pred_ts = timeseries(modified_ClearSkyGHI_ineichen,datestr(dt_next_window));
figure(1);
subplot(3,1,1)
hold on;
plot(CS_pred_ts,':k')


%%

% 
% 
% 
% 
% 
% 
% 
% %------------- future windows
% 
% % Variables have different time stamps
% % Resampling variables / to make them all uniform
% 
% now_vec = datevec(now); 
% %now_vec = datevec(dt_GHI(end))
% %datestr(now_vec)
% 
% day_now = datenum(now_vec(1:3));
% %datestr(day_now) % debug
% hr_now  =  datenum(now_vec(4));
% min_now = datenum(now_vec(5)) ;
% 
% % 10 min interval 1min resolution
% next_ten_min  = 1 /24/60:1 /24/60:10 /24/60 ;
% 
% % Up to minuntes
% now_time = day_now + hr_now/24 + min_now/24/60; 
% %datestr(now_time) % Debug
% 
% % Past window
% dt_future_window = now_time + next_ten_min';
% %datestr(dt_past_window) % Debug
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % -------------- Apply Clear Sky Model for the Future --------------------------
% Location.latitude =  8+16/60+38.50/60/60;     %   8�16'38.50" E   Local Longitude of CHCRC.C1 roof (negative if west)
% Location.longitude = 47+27/60+33.92/60/60;    %   47�27'33.92" N  Local latitude of CHCRC.C1 roof  (positive if north of equator)
% Location.altitude =  400; %455;                     %   [m]  
% 
% 
%              
% DN = dt_future_window; %GHIData{iter}.Time %+  start_date    ;
% 
% % Siwtzerland has a Central European Summer Time (CEST) +0200 UTC 
% UTC = 2; %4.5; %  UTC is 2 but I have noticed an offset of 2:30 min.
% Time = pvl_maketimestruct(DN, UTC);
% 
% % Clear sky GHI, DNI, and DHI from Ineichen/Perez model
% [ClearSkyGHI_ineichen, ClearSkyDNI_ineichen, ClearSkyDHI_ineichen]= pvl_clearsky_ineichen(Time, Location);       
% 
% cloudyness_factor = 0.8 * mean(ClearSkyGHI_ineichen ./ GHI_past)  ;
% modified_ClearSkyGHI_ineichen= ClearSkyGHI_ineichen ./ cloudyness_factor;
% 
% % Clear sky GHI from Haurwitz model
% [SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(Time, Location);
% ApparentZenith = 90-ApparentSunEl;
% 
% ClearSkyGHI_haurwitz = pvl_clearsky_haurwitz(ApparentZenith);
% 
% 
% 
% 
% 
% 
% 
% 
% figure(2) % Use the same figure
% 
% subplot(3,2,2)
% hold on, plot(dt_future_window,modified_ClearSkyGHI_ineichen,'r'),hold off
% datetick('x','HH:MM','keeplimits')
% xlim([dt_future_window(1) - 1/24/60, dt_future_window(end) + 1/24/60  ])
% title(['10 min Future Window of GHI for ', datestr(floor(dt_GHI(end)))])
% ylabel('W/m^2')
% 
% 
% % Temperature
% %subplot(3,2,4)
% % plot(dt_past_window,Temp_past),
% % xlim([dt_past_window(1) - 1/24/60, dt_past_window(end) + 1/24/60  ])
% % ylabel ('�C')
% % title(['10 min Past Window of Temperature for ',datestr(floor(dt_GHI(end)))])
% % datetick('x','HH:MM','keeplimits')
% 
% %Power
% %subplot(3,1,6)
% % plot(dt_past_window,P_PV_past), 
% % xlim([dt_past_window(1) - 1/24/60, dt_past_window(end) + 1/24/60  ])
% % ylabel ('W')
% % title(['10 min Past Window of Power for ', datestr(floor(dt_PV(end)))])
% % datetick('x','HH:MM','keeplimits')
% % %set(gca,'XTick',0:0.2:1); %Set the XTicks values
% % %datetick('x','HH','keepticks')
% % %xlim([0 1])
% % %ylim([0 1.2e3])
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
