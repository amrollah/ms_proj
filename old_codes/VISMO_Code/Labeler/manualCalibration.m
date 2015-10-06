%% put the file that you want to use and 
im_dir = './Calibration_Sets/3/';
LabelerGUI(LabelerPrep('dir',im_dir));

%%
% after the gui is shot down, go and load the labelers and run this part...
%load 'calib_modelGood.mat';
load 'calib_model_01_2014.mat';
colors_num_class = {1,2,3,4,5,6,7,8,9};
colors_num_keys = {'m', 'g', 'y', 'b', 'p', 'r', 'k', 'w', 'c'};
clr_num_map = containers.Map(colors_num_keys(:),colors_num_class(:));
label_dataset = [];
CH_Chim_Theta = 270+(180-154.1431);   % 154+(08/60)+(35/3600) at this bearing wrt y axis it is found, converted to bearing wrt x axis for polar coordinate use...
ABB_Chim_Theta = 270+(180-186.5669);  % 186+(34/60)+(01/3600) at this bearing wrt y axis it is found, converted to bearing wrt x axis for polar coordinate use...
correction_theta = [];
num_samples = 0;
for i = length(im_labelers):-1:1
    for j = 1:length(im_labelers(i).im_labels)
        num_samples = num_samples + 1;
        l_str = im_labelers(i).im_labels(j);
        %box_pix_pose = l_str.bbox(2:-1:1)+(l_str.bbox(4:-1:3)/2);
        box_pix_pose = l_str.bbox(1:2)+(l_str.bbox(3:4)/2);
        uni_sph_pose = cam2world(box_pix_pose(2:-1:1)',ocam_model);
        theta1 = cart2pol(-uni_sph_pose(1),uni_sph_pose(2))*(180/pi);
        if(clr_num_map(l_str.class)==5) % this is the Hospital's thing (dark blue)
            theta2 = CH_Chim_Theta;
            correction_theta(num_samples) = mod(180-(theta1-theta2),360); 
        elseif(clr_num_map(l_str.class)==3) % this is the ABB's thing (yellow)
            theta2 = ABB_Chim_Theta;
            correction_theta(num_samples) = mod(180-(theta1-theta2),360);
        elseif(clr_num_map(l_str.class)==7) % this is the construction thing
        end
    end
end

theta = mean(correction_theta);
ecam_model.theta = theta;
ecam_model.R = [cos(theta*pi/180),cos((theta-90)*pi/180);...
                sin(theta*pi/180),sin((theta-90)*pi/180)];
ecam_model.Rinv = inv([cos(theta*pi/180),cos((theta-90)*pi/180);...
                sin(theta*pi/180),sin((theta-90)*pi/180)]);
cd ../Utils/Data
save('extrinsic_calib_model6', 'ecam_model');
cd ../../Labeler