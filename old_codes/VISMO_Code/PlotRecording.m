%% Prepared for Andi for him to be able to see the performance of the 
% tracking algorithm and the kalman filter. 
% B. Zeydan, 22. May. 2013

%%
e_val = exist('recordings');
if ~(e_val == 7)
    error('no recording exists');
end
cd recordings;
c_list = sort(cell2mat(cellfun(@(x) str2num(x),cellstr(ls),'UniformOutput',false)));

if ~isempty(c_list)
    b_num = c_list(end);
else
    error('no recording exists');
end

b_str = num2str(b_num);
cd(b_str);
e_val1 = exist('recording_obs');
e_val2 = exist('recording_pred');
e_val3 = exist('db_objs');
e_val4 = exist('db_ass');
if((e_val1 == 0) || (e_val2 == 0) || (e_val3 == 0) || (e_val4 == 0))
    e_val1 = exist('observed.mat');
    e_val2 = exist('predicted.mat');
    e_val3 = exist('observations.mat');
    e_val4 = exist('assignments.mat');
    if((e_val1 ~= 2) || (e_val2 ~= 2) || (e_val3 ~= 2) || (e_val4 ~= 2))
        error('incomplete recording please delete the latest unsuccessful recording');
    else
        load('observed.mat');
        load('predicted.mat');
        load('observations.mat');
        load('assignments.mat');
    end
else
    warning('Already existing recordings. Not loading. To reload delete recording_obs, recording_preds variables from your workspace and try again.');
end
cd ../..

%%
for i = 1:length(recording_obs)
    temp_ob = recording_obs{i};
    if(length(temp_ob)~=0)
        centroid_obs(i,[temp_ob(:).id],:) = cell2mat(arrayfun(@(y) y.centroid ,temp_ob(:),'UniformOutput',false));
        centroid_obs_proj(i,[temp_ob(:).id],:) = cell2mat(arrayfun(@(y) y.centroidProj ,temp_ob(:),'UniformOutput',false));
        speed_obs(i,[temp_ob(:).id],:) = cell2mat(arrayfun(@(y) y.obsSpeed ,temp_ob(:),'UniformOutput',false));
    end
end

for i = 1:length(recording_pred)
    temp_pr = recording_pred{i};
    if(length(temp_pr)~=0)
        centroid_pred(i,[temp_pr(:).id],:) = cell2mat(arrayfun(@(y) y.predCentroid ,temp_pr(:),'UniformOutput',false));
        centroid_pred_proj(i,[temp_pr(:).id],:) = cell2mat(arrayfun(@(y) y.predCentroidProj' ,temp_pr(:),'UniformOutput',false));
        speed_pred(i,[temp_pr(:).id],:) = cell2mat(arrayfun(@(y) y.predSpeed ,temp_pr(:),'UniformOutput',false));
    end 
end
%cellfun(@(x) arrayfun(@(y) y.id, x,'UniformOutput',false),recording_obs(:),'UniformOutput',false);


%%

i = 40; % the number of the tracks should be less than the total number of tracks

if(i<=size(centroid_pred,2))
    figure; subplot(2,2,1);
    plot(1:length(centroid_obs(:,i,1)),centroid_obs(:,i,1),'r');
    hold on;
    plot(1:length(centroid_pred(:,i,1)),centroid_pred(:,i,1),'k');
    title('X prediction in black, observation in red')
    subplot(2,2,2);  plot(1:length(centroid_obs(:,i,2)),centroid_obs(:,i,2));
    hold on;
    plot(1:length(centroid_pred(:,i,2)),centroid_pred(:,i,2),'m');
    title('Y prediction in magenta, observation in blue')
    subplot(2,2,3);
    plot(1:length(speed_obs(:,i,1)),speed_obs(:,i,1),'r');
    hold on;
    plot(1:length(speed_pred(:,i,1)),speed_pred(:,i,1),'k');
    title('Xdot prediction in black, observation in red')
    subplot(2,2,4);  plot(1:length(speed_obs(:,i,2)),speed_obs(:,i,2));
    hold on;
    plot(1:length(speed_pred(:,i,2)),speed_pred(:,i,2),'m');
    title('Ydot prediction in magenta, observation in blue');
else
    warning(['No Such Recording. Nothing to plot. The i should be between 1 and '  num2str(size(centroid_pred,2))]);
end