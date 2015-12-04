%close all; 
clear all; 
%clc;

%% settings
addpath(genpath('C:\Users\chamsei\Documents\GitHub\ms_proj\old_codes\vismolib'))
addpath(genpath('C:\Users\chamsei\Documents\GitHub\ms_proj\old_codes\solvers'));
vmlconfig_cavriglia;
conf = evalin('base','VMLCONF');
save_figs = true;
show_power = true;
shift_power = false;
text_labels = false;

%% Path of data
% folder of McClear export files
mcclear_path = 'C:\data\mcclear\';
% base folder for Cavriglia image and data files
cav_data_path = 'U:\HDRImages_Cavriglia\img\';
local_cav_data_path = 'C:\data\cav\log_data\';

%% find all available days for irradiation visualization
days = dir(cav_data_path); % get all content of Cavriglia image folder
days = days(vertcat(days.isdir)); % filter only folders
days = days(8:end-1)'; % skip first 5 folders (first two are . and .., data of next 5 are corrupted) and skip last folder, because log data is not still transfered for current day.
%'2015_07_19'

invalid_dates = {};
counter = 1;
%%% loop over all days
for d=1:size(days,2)
    date = days(d).name;
    date = '2015_08_03';
    obj = vmlSeq(date,[6 20]);
    
    disp(date);
    if ~isempty(obj.P)
        P_tmp = sort(obj.P(:,2));
    end
    if isempty(obj.P) || mean(P_tmp(1:500))< -1000
        disp('Invalid power log files.');
        invalid_dates{counter} = date;
        counter = counter + 1;
    elseif show_power
        Irr_ts = timeseries(obj.Irr(:,2:end), obj.Irr(:,1));
        Temp_ts = timeseries(obj.Temp(:,2), obj.Temp(:,1));
        Irr = resample(Irr_ts, obj.P(:,1)); 
        temp = resample(Temp_ts, obj.P(:,1));
        pm=floor(size(obj.P,1)/2);
        
        txt_idx_x = 1:50:(numel(Irr.data(:,1))-4);
        txt_idx_y = 1:50:(numel(obj.P(:,2))-4);
        
        %% 
        if shift_power
            min_shift = -5;
            max_shift = 30;
            figure;
            sbx = floor(sqrt((max_shift-min_shift)/2));
            sby = ceil((max_shift-min_shift)/(2*sbx));
            sh_counter = 1;
            for shift=min_shift:2:max_shift-1
                subplot(sbx,sby,sh_counter);
                if shift < 0
                    plot(Irr.data((1-shift):end,1), obj.P(1:(end+shift),2),'b.');
                else
                    plot(Irr.data(1:(end-shift),1), obj.P((1+shift):end,2),'b.');
                end
                title(strcat(num2str(shift), ' shift'));
                sh_counter=sh_counter+1;
            end
        end        
        
        %% regreesion
        lin1_idxs = find(Irr.data(:,1)>=140);
        lin2_idxs = find(Irr.data(:,1)<140);
        X1 = [ones(size(Irr.data(lin1_idxs,1))), Irr.data(lin1_idxs,1), temp.data(lin1_idxs,1), Irr.data(lin1_idxs,1).*temp.data(lin1_idxs,1)];
%         b1 = regress(obj.P(lin1_idxs,2),X1);
        b1 = fitlm(X1,obj.P(lin1_idxs,2),'RobustOpts','on');
        
        X2 = [ones(size(Irr.data(lin2_idxs,1))), Irr.data(lin2_idxs,1), Irr.data(lin2_idxs,1).^2, temp.data(lin2_idxs,1), Irr.data(lin2_idxs,1).*temp.data(lin2_idxs,1)];
%         b2 = regress(obj.P(lin2_idxs,2),X2);
        b2 = fitlm(X2,obj.P(lin2_idxs,2),'RobustOpts','on');
        
        est_pw = zeros(size(obj.P,1),2);
        est_pw(:,1) = obj.P(:,1);
        est_pw(lin1_idxs,2) = b1.Fitted;%X1*b1;
        est_pw(lin2_idxs,2) = b2.Fitted;%X2*b2;
        
        %%
        figure;
        subplot(1,2,1);
        plot(Irr.data(1:end-pm,1), obj.P(1:end-pm,2),'b.');
        hold on;
        plot(Irr.data(end-pm:end,1), obj.P(end-pm:end,2),'r.');
        hold on;
        plot(Irr.data(:,1), est_pw(:,2),'k.');
        
        title(strrep(date, '_', '/'));
        xlabel('Irradiance');
        ylabel('Power');
        if text_labels
            [~,~,~,labels(:,1),labels(:,2),labels(:,3)] = datevec(obj.P(txt_idx_x,1));
            labels(:,3) = floor(labels(:,3));
            h = text(Irr.data(txt_idx_x,1),obj.P(txt_idx_y,2),strcat(num2str(labels(:,1)),':',num2str(labels(:,2)),':',num2str(labels(:,3))),'FontSize',5);
            set (h, 'Clipping', 'on');
        end
        hold off;
        
        subplot(1,2,2);
        xtick = 1:100:size(Irr.time);
        [~,~,~,x_lables(:,1),x_lables(:,2),x_lables(:,3)] = datevec(Irr.time(xtick));
        x_lables(:,3) = floor(x_lables(:,3));
        
        plot(Irr.data(:,1), '-g');
        hold on;
%         plot(obj.ClearSkyRef(:,2), '--b');
%         hold on;
%         plot(obj.ClearSkyOrigRef(:,2), '--r');
       
        est_power = linmap(est_pw(:,2),[min(Irr.data(:,1)),max(Irr.data(:,1))]);
        plot(est_power, '-r');
        hold on;
        t_power = linmap(obj.P(:,2),[min(Irr.data(:,1)),max(Irr.data(:,1))]);
        plot(t_power, '-k');
        hold on;
        plot(10*temp.data(:,1), '-c');
        
        set(gca, 'XTick', xtick, 'XTickLabel',strcat(num2str(x_lables(:,1)),':',num2str(x_lables(:,2)),':',num2str(x_lables(:,3))));
        rotateXLabels(gca, 90);
        title(strcat(strrep(date, '_', '/'), ' irradiance'));
        xlabel('Time');
        ylabel('Irradiance and Power');
        
        pause(1);
    end
end
