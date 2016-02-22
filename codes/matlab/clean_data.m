clear all;
close all;
clc;
img_save_path='';
prj_path='';
proj_path; 
load('calc\img_data4.mat', 'data');
cl_data = {};
for i=1:length(data)
    d = data{i};
    img = imread([img_save_path d.day '__' num2str(d.j) '.jpeg']);
    red_m=mean(mean(img(:,:,1)));
    if red_m>100
       disp([num2str(i), '     ', num2str(red_m)]);
    end
    if red_m<150 && d.diff_irr > .05*d.irr(1) && d.diff_irr < .6*d.clear_irr(2)
%         figure(1);
%         imshow(img);
        cl_data{end+1}=d;
    end
end
 
save('calc\data_clean.mat', 'cl_data');