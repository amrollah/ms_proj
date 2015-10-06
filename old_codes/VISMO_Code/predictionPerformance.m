
cc = 75;
algIrrPerfs = zeros(length(perfs),cc);
persIrrPerfs2 = zeros(length(perfs),cc);
persIrrPerfs = zeros(length(perfs),cc);
csIrrPerfs = zeros(length(perfs),cc);
algOccPerfs = zeros(length(perfs),cc);
persOccPerfs = zeros(length(perfs),cc);
cover100 = zeros(length(perfs),cc);
cover500 = zeros(length(perfs),cc);
cover1000 = zeros(length(perfs),cc);
% cover100 = zeros(length(perfs),1);
% cover500 = zeros(length(perfs),1);
% cover1000 = zeros(length(perfs),1);
for i = 1:length(perfs)
    algIrrPerfs(i,:) = abs(perfs(i).mesIrr - perfs(i).algIrr);
    persIrrPerfs(i,:) = abs(perfs(i).mesIrr - perfs(i).persIrr);
    csIrrPerfs(i,:) = abs(perfs(i).mesIrr - perfs(i).csIrr);
    %algOccPerfs(i,:) = ~xor(perfs(i).mesOcc, perfs(i).algOcc);
    %persOccPerfs(i,:) = ~xor(perfs(i).mesOcc, perfs(i).persOcc); 
    algOccPerfs(i,:) = 1-abs((perfs(i).mesOcc - perfs(i).algOcc));
    persOccPerfs(i,:) = 1-abs(perfs(i).mesOcc - perfs(i).persOcc);
    if(i>1)
        persIrrPerfs2(i,:) = abs(perfs(i).mesIrr - perfs(i-1).mesIrr(1));
    end
    cover100(i,:) = perfs2(i).mesCov100;
    cover500(i,:) = perfs2(i).mesCov500;
    cover1000(i,:) = perfs2(i).mesCov1000;
    cover101(i,:) = perfs(i).mesCov100;
    cover501(i,:) = perfs(i).mesCov500;
    cover1001(i,:) = perfs(i).mesCov1000;
end
occs = cat(1,perfs.mesOcc);
mean(occs(:,2))
%algIrrPerfs(:,round(2.5:2.5:end))
%
plot([0.5:0.5:15],mean(persIrrPerfs(:,round(2.5:2.5:end))),[0.5:0.5:15],mean(algIrrPerfs(:,round(2.5:2.5:end))));
figure, plot([0.5:0.5:15],sqrt(mean(persIrrPerfs2(:,round(2.5:2.5:end)).^2)),[0.5:0.5:15],sqrt(mean(algIrrPerfs(:,round(2.5:2.5:end)).^2)));
plot(1:length(perfs),persIrrPerfs(:,end),1:length(perfs),algIrrPerfs(:,end));
%plot(1:length(perfs),persOccPerfs(:,end-15),1:length(perfs),algOccPerfs(:,end-15))
% mean(persIrrPerfs(:,round(2.5:2.5:end)))
% median(persIrrPerfs(:,round(2.5:2.5:end)))
% min(persIrrPerfs(:,round(2.5:2.5:end)))
% 
% mean(algIrrPerfs(:,round(2.5:2.5:end)))
% median(algIrrPerfs(:,round(2.5:2.5:end)))
% min(algIrrPerfs(:,round(2.5:2.5:end)))
% sqrt(mean(algIrrPerfs(:,round(2.5:2.5:end)).^2))
% sqrt(mean(persIrrPerfs(:,round(2.5:2.5:end)).^2))
% 
% mean(csIrrPerfs(:,round(2.5:2.5:end)))
% 
% mean(algOccPerfs(:,round(2.5:2.5:end)))
% mean(persOccPerfs(:,round(2.5:2.5:end)))

mean(persIrrPerfs(:,:))
median(persIrrPerfs(:,:))
min(persIrrPerfs(:,:))

mean(persIrrPerfs2(2:end,:))
median(persIrrPerfs2(2:end,:))
min(persIrrPerfs2(2:end,:))

mean(algIrrPerfs(:,:))
median(algIrrPerfs(:,:))
min(algIrrPerfs(:,:))
mean(algIrrPerfs(:,round(7.5:7.5:end)))

sqrt(mean(algIrrPerfs(:,:).^2))
sqrt(mean(persIrrPerfs(:,:).^2))
sqrt(mean(persIrrPerfs2(2:end,:).^2))
sqrt(mean(csIrrPerfs(:,:).^2))

mean(csIrrPerfs(:,:))

mean(algOccPerfs(:,:))
mean(persOccPerfs(:,:))

mesIrr = cat(1,perfs.mesIrr);
algIrr = cat(1,perfs.algIrr);
persIrr = cat(1,perfs.persIrr);
csIrr = cat(1,perfs.csIrr);
mesOcc = cat(1,perfs.mesOcc);
algOcc = cat(1,perfs.algOcc);
time_s = cat(1,perfs.frameDateNum);

% tsS5 = timeseries(sensor5Data(:,14),datestr(sensor5Data(:,1)));
% tsS5 = setinterpmethod(tsS5,'zoh');
% tsS5_n = resample(tsS5,datestr(time_s));

figure, plot(time_s,(csIrr(:,1).*(1-cover101))/8,time_s,(csIrr(:,1).*(1-cover501))/8,time_s,(csIrr(:,1).*(1-cover1001))/8, 'LineWidth',1.5);
set(gca,'XTick',[time_s(1):2/(60*24):time_s(end)]);
datetick('x','keepticks');
set(gca,'FontSize',14);
legend('100x100', '500x500','1000x1000','Location','Southwest');

predHorizon = 25;
timeHorizon = predHorizon*4*3/(60*60*24);
figure, plot(time_s+timeHorizon,(csIrr(:,predHorizon).*(1-cover100(:,predHorizon)))/8,time_s+timeHorizon,(csIrr(:,predHorizon).*(1-cover500(:,predHorizon)))/8,time_s+timeHorizon,(csIrr(:,predHorizon).*(1-cover1000(:,predHorizon)))/8, 'LineWidth',1.5);
set(gca,'XTick',[time_s(1)+timeHorizon:2/(60*24):time_s(end)+timeHorizon]);
datetick('x','keepticks');
set(gca,'FontSize',14);
legend('100x100', '500x500','1000x1000','Location','Southwest');

predHorizon = 75;
timeHorizon = predHorizon*4*3/(60*60*24);
figure, plot(time_s,(csIrr(:,1).*(1-cover101))/8,time_s+timeHorizon/3,(csIrr(:,predHorizon/3).*(1-cover100(:,predHorizon/3)))/8,time_s+timeHorizon*(2/3),(csIrr(:,predHorizon*(2/3)).*(1-cover100(:,predHorizon*(2/3))))/8,time_s+timeHorizon,(csIrr(:,predHorizon).*(1-cover100(:,predHorizon)))/8, 'LineWidth',1.5);
set(gca,'XTick',[time_s(1):2/(60*24):time_s(end)+timeHorizon]);
datetick('x','keepticks');
set(gca,'FontSize',14);
legend('100x100 current', '100x100 5 mins','100x100 10 mins','100x100 15 mins','Location','Southwest');

figure, plot(time_s,(csIrr(:,1).*(1-cover501))/8,time_s+timeHorizon/3,(csIrr(:,predHorizon/3).*(1-cover500(:,predHorizon/3)))/8,time_s+timeHorizon*(2/3),(csIrr(:,predHorizon*(2/3)).*(1-cover500(:,predHorizon*(2/3))))/8,time_s+timeHorizon,(csIrr(:,predHorizon).*(1-cover500(:,predHorizon)))/8, 'LineWidth',1.5);
set(gca,'XTick',[time_s(1):2/(60*24):time_s(end)+timeHorizon]);
datetick('x','keepticks');
set(gca,'FontSize',14);
legend('500x500 current', '500x500 5 mins','500x500 10 mins','500x500 15 mins','Location','Southwest');

timeHorizon = 0;
figure, plot(time_s,(csIrr(:,1).*(1-cover1001))/8,time_s+timeHorizon/3,(csIrr(:,1).*(1-cover1000(:,predHorizon/3)))/8,time_s+timeHorizon*(2/3),(csIrr(:,1).*(1-cover1000(:,predHorizon*(2/3))))/8,time_s+timeHorizon,(csIrr(:,1).*(1-cover1000(:,predHorizon)))/8, 'LineWidth',1.5);
set(gca,'XTick',[time_s(1):2/(60*24):time_s(end)+timeHorizon]);
datetick('x','keepticks');
set(gca,'FontSize',14);
legend('1000x1000 current', '1000x1000 5 mins','1000x1000 10 mins','1000x1000 15 mins','Location','Southwest');

%figure, plot(1:length(perfs),mesIrr(:,25),1:length(perfs),algIrr(:,25));
figure, plot(time_s, mesIrr(:,25), time_s, algIrr(:,25), time_s, csIrr(:,25), 'LineWidth',1.5);
set(gca,'XTick',[time_s(1):1/12:time_s(end)]);
datetick('x',15);
set(gca,'FontSize',14);
legend('Measured', 'Predicted','Clear Sky','Location','Southwest');
xlabel('Time');
ylabel('Irradiance in (W/m^2)');

figure, plot(time_s, mesOcc(:,25), time_s, algOcc(:,25), 'LineWidth',1.5);
set(gca,'XTick',[time_s(1):1/12:time_s(end)]);
datetick('x',15);
set(gca,'FontSize',14);
legend('Measured', 'Predicted','Location','Southwest');
xlabel('Time');
ylabel('Binary Occlusion');

figure, plot(time_s, mesIrr(:,end), 'LineWidth',1.5);
set(gca,'XTick',[time_s(1):1/12:time_s(end)]);
datetick('x',15);
set(gca,'FontSize',14);`
%legend('Measured', 'Predicted','Clear Sky','Location','Southwest');
xlabel('Time');
ylabel('Irradiance in (W/m^2)');


%%
%perfs_st = perfs;
%perfs = perfs_st(1:end/2);
cc = 25;

load 'PVGHIData_plus_CSK.mat';
selected_month = 8;
selected_year =  2013;
selected_day = 12;

idx=[];

for i=1:numel(GHIData)
     [y, m, d] = datevec (GHIData{i}.DateNumber);
      if selected_month == m && selected_year == y && selected_day == d      
          idx = i;      
      end     
end
d=datenum([2013,11,26,12,0,0])-GHIData{idx}.DateNumber;
offset = datenum([0,0,0,0,33,2])-datenum([0,0,0,0,0,d*12.3]);
timeSeriesPV = timeseries(GHIData{idx}.PV_output(:),datestr(GHIData{idx}.Time+GHIData{idx}.DateNumber+offset));
timeSeriesTemp = timeseries(GHIData{idx}.ambTemp(:),datestr(GHIData{idx}.Time+GHIData{idx}.DateNumber+offset));
powerDB = resample(timeSeriesPV,datestr([perfs_st(:).frameDateNum]'));
tempDB = resample(timeSeriesTemp,datestr([perfs_st(:).frameDateNum]'));

sampleData = struct(...
    'dateNum', {},...
    'predIrr', {},...
    'measIrr', {},...
    'measPV', {},...
    'measTemp', {});

for i = 1:length(perfs)
    newData = struct(...
        'dateNum', perfs(i).frameDateNum,...
        'predIrr', perfs(i).algIrr,...
        'measIrr', perfs(i).mesIrr,...
        'measPV', powerDB.Data(i+1:i+25)',...
        'measTemp', tempDB.Data(i+1:i+25)');
    sampleData(end+1) = newData;
end