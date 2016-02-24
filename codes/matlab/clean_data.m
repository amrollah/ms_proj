clear all;
close all;
clc;
img_save_path='';
prj_path='';
proj_path; 
load('calc\data_clean.mat', 'data');
load('calc\max_irr.mat', 'values');
cl_data = {};
for i=1:length(data)
    d = data{i};
%     img = imread([img_save_path d.day '__' num2str(d.j) '.jpeg']);
%     red_m=mean(mean(img(:,:,1)));
%     if red_m>100
%        disp([num2str(i), '     ', num2str(red_m)]);
%     end
%     if red_m<150 && d.diff_irr > .05*d.irr(1) && d.diff_irr < .6*d.clear_irr(2)
%         figure(1);
%         imshow(img);
%     if ~(strcmp(d.day,'2015_09_15') || strcmp(d.day,'2015_09_16') || strcmp(d.day,'2015_09_17') || strcmp(d.day,'2015_09_18')...
%             || strcmp(d.day,'2015_10_27') || strcmp(d.day,'2015_10_28') || strcmp(d.day,'2015_10_29') || strcmp(d.day,'2015_10_30')...
%             || strcmp(d.day,'2015_10_31') || strcmp(d.day,'2015_11_01') || strcmp(d.day,'2015_11_02') || strcmp(d.day,'2015_11_03')...
%             || strcmp(d.day,'2015_11_04') || strcmp(d.day,'2015_11_05') || strcmp(d.day,'2015_11_06') || strcmp(d.day,'2015_11_07')...
%             || strcmp(d.day,'2015_11_08') || strcmp(d.day,'2015_11_09') || strcmp(d.day,'2016_01_02'))
%         cl_data{end+1}=d;
%     end
max_irr = NaN;
for j=1:length(values)
    if strcmp(values{j}.date, date)
        max_irr= values{j}.max_irr;
    end
end

if d.clear_irr(1)/max_irr < .1
   disp(datevec(d.time)); 
end
end
% data=cl_data;
% save('calc\data_clean2.mat', 'data');