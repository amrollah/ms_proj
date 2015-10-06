%% load the necessary objects
version2 = 1;
version3 = 0;
current_folder = pwd;
cd(save_folder);
load('details.mat');
load('labels.mat');
cd(current_folder);
if(version2 || version3)
    load 'calib_model_01_2014.mat';
    load 'extrinsic_calib_model_01_2014.mat';
    load 'model3D_Bigger.mat';
    [X,Y] = meshgrid(1:ocam_model.height,1:ocam_model.width);
    deg2radf = pi/180;
end
e_val = exist('im_arr');
if(e_val ~= 1)
    im_arr = cellfun(@(x) imread([im_str.im_dir x]),im_str.im_list,'UniformOutput',false);
    name_arr = cellfun(@(x) strtok(char(x),'.'),im_str.im_list,'UniformOutput',false);
    date_arr = cell2mat(cellfun(@(x) strread(x,'%d','delimiter','_'),name_arr,'UniformOutput',false)')';
    date_arr = date_arr(:,1:6);
end
%%
colors_num_class = {1,2,3,4,5,6,7,8,9};
colors_num_keys = {'m', 'g', 'y', 'b', 'p', 'r', 'k', 'w', 'c'};
clr_num_map = containers.Map(colors_num_keys(:),colors_num_class(:));
label_dataset = [];
% im_arr2 = im_arr;
if(version2 || version3)
    for k = 1:length(im_arr)
        im_c = cell2mat(im_arr(k));
        imDateVec = date_arr(k,:);
        DN = datenum(imDateVec);
        Time = pvl_maketimestruct(DN, model3D.UTC);
        [SunAz, ~, ApparentSunEl, ~]=pvl_ephemeris(Time, model3D.Location);
        zenith = 90-ApparentSunEl;
        azimuth = SunAz;
        [x,y,z] = sph2cart((90-azimuth)*deg2radf,(90-zenith)*deg2radf,1);
        sunPositionReal = [x;y;-z];
        sunPosition = world2cam(inv(ecam_model.R)*sunPositionReal,ocam_model);
        sunPatch = sqrt((X-sunPosition(2)).^2+(Y-sunPosition(1)).^2)<=180;
        r_b = double(im_c(:,:,1))./double(im_c(:,:,3)); % R/B thresholding applied
        rgb = double(im_c(:,:,1))+double(im_c(:,:,2))+double(im_c(:,:,3)); % to distinguish other objects like buildings from clouds, a test on brightness
        mask = r_b > 0.72 & rgb > 50;             %% brightness auto thing is needed because of this ..
        mask_arr{k} = mask;
        g_b = double(im_c(:,:,2))./double(im_c(:,:,3));
        mask2 = g_b > 0.77;
        mask(sunPatch) = mask2(sunPatch);
        mask_arr2{k} = mask;
        sun_pos(k,:) = sunPosition';
    end
end
for i = length(im_labelers):-1:1
    for j = 1:length(im_labelers(i).im_labels)
        l_str = im_labelers(i).im_labels(j);
        im_c = cell2mat(im_arr(l_str.id));
        im_c1 = double(im_c(:,:,1));
        im_c2 = double(im_c(:,:,2));
        im_c3 = double(im_c(:,:,3));
        if(version2)
            mask = cell2mat(mask_arr2(l_str.id));
            % 1  2  3   4    5    6   7   8   9   10     11       12
            % r, g, b, r_b, r_g, g_b r-b r-g g-b r+g+b  label  decision
            if(size(label_dataset,1) == 0)
                label_dataset(1:sum(sum(l_str.mask)),:) = [im_c1(l_str.mask),im_c2(l_str.mask),im_c3(l_str.mask), ...
                    im_c1(l_str.mask)./im_c3(l_str.mask),im_c1(l_str.mask)./im_c2(l_str.mask),im_c2(l_str.mask)./im_c3(l_str.mask),...
                    im_c1(l_str.mask)-im_c3(l_str.mask),im_c1(l_str.mask)-im_c2(l_str.mask),im_c2(l_str.mask)-im_c3(l_str.mask),...
                    im_c1(l_str.mask)+im_c2(l_str.mask)+im_c3(l_str.mask),repmat(clr_num_map(l_str.class),sum(sum(l_str.mask)),1),mask(l_str.mask)];
            else
                label_dataset(end:end+sum(sum(l_str.mask))-1,:) = [im_c1(l_str.mask),im_c2(l_str.mask),im_c3(l_str.mask), ...
                    im_c1(l_str.mask)./im_c3(l_str.mask),im_c1(l_str.mask)./im_c2(l_str.mask),im_c2(l_str.mask)./im_c3(l_str.mask),...
                    im_c1(l_str.mask)-im_c3(l_str.mask),im_c1(l_str.mask)-im_c2(l_str.mask),im_c2(l_str.mask)-im_c3(l_str.mask),...
                    im_c1(l_str.mask)+im_c2(l_str.mask)+im_c3(l_str.mask),repmat(clr_num_map(l_str.class),sum(sum(l_str.mask)),1),mask(l_str.mask)];
            end
%             im_arr2{l_str.id} = objectAnnotation(double(cell2mat(im_arr2(l_str.id))),'Rectangle',round(l_str.bbox),l_str.class,10);
        elseif(version3)
            mask = cell2mat(mask_arr2(l_str.id));
            % 1  2  3   4    5    6   7   8   9   10       11       12       13
            % r, g, b, r_b, r_g, g_b r-b r-g g-b r+g+b dist_to_sun label  decision
            [x_pos,y_pos] = ind2sub(size(l_str.mask),find(l_str.mask));
            sun_dists = sqrt(((x_pos-sun_pos(l_str.id,1))./size(l_str.mask,1)).^2 + ((y_pos-sun_pos(l_str.id,2))./size(l_str.mask,2)).^2);
            if(size(label_dataset,1) == 0)
                label_dataset(1:sum(sum(l_str.mask)),:) = [im_c1(l_str.mask),im_c2(l_str.mask),im_c3(l_str.mask), ...
                    im_c1(l_str.mask)./im_c3(l_str.mask),im_c1(l_str.mask)./im_c2(l_str.mask),im_c2(l_str.mask)./im_c3(l_str.mask),...
                    (im_c1(l_str.mask)-im_c3(l_str.mask))/255,(im_c1(l_str.mask)-im_c2(l_str.mask))/255,(im_c2(l_str.mask)-im_c3(l_str.mask))/255,...
                    (im_c1(l_str.mask)+im_c2(l_str.mask)+im_c3(l_str.mask))/765,sun_dists,repmat(clr_num_map(l_str.class),sum(sum(l_str.mask)),1),mask(l_str.mask)];
            else
                label_dataset(end:end+sum(sum(l_str.mask))-1,:) = [im_c1(l_str.mask),im_c2(l_str.mask),im_c3(l_str.mask), ...
                    im_c1(l_str.mask)./im_c3(l_str.mask),im_c1(l_str.mask)./im_c2(l_str.mask),im_c2(l_str.mask)./im_c3(l_str.mask),...
                    (im_c1(l_str.mask)-im_c3(l_str.mask))/255,(im_c1(l_str.mask)-im_c2(l_str.mask))/255,(im_c2(l_str.mask)-im_c3(l_str.mask))/255,...
                    (im_c1(l_str.mask)+im_c2(l_str.mask)+im_c3(l_str.mask))/765,sun_dists,repmat(clr_num_map(l_str.class),sum(sum(l_str.mask)),1),mask(l_str.mask)];
            end
        else
            %im_c = imread(im_str.im_list(l_str.id));
            % 1  2  3   4    5    6   7   8   9   10    11
            % r, g, b, r_b, r_g, g_b r-b r-g g-b r+g+b label
            if(size(label_dataset,1) == 0)
                label_dataset(1:sum(sum(l_str.mask)),:) = [im_c1(l_str.mask),im_c2(l_str.mask),im_c3(l_str.mask), ...
                    im_c1(l_str.mask)./im_c3(l_str.mask),im_c1(l_str.mask)./im_c2(l_str.mask),im_c2(l_str.mask)./im_c3(l_str.mask),...
                    im_c1(l_str.mask)-im_c3(l_str.mask),im_c1(l_str.mask)-im_c2(l_str.mask),im_c2(l_str.mask)-im_c3(l_str.mask),...
                    im_c1(l_str.mask)+im_c2(l_str.mask)+im_c3(l_str.mask),repmat(clr_num_map(l_str.class),sum(sum(l_str.mask)),1)];
            else
                label_dataset(end:end+sum(sum(l_str.mask))-1,:) = [im_c1(l_str.mask),im_c2(l_str.mask),im_c3(l_str.mask), ...
                    im_c1(l_str.mask)./im_c3(l_str.mask),im_c1(l_str.mask)./im_c2(l_str.mask),im_c2(l_str.mask)./im_c3(l_str.mask),...
                    im_c1(l_str.mask)-im_c3(l_str.mask),im_c1(l_str.mask)-im_c2(l_str.mask),im_c2(l_str.mask)-im_c3(l_str.mask),...
                    im_c1(l_str.mask)+im_c2(l_str.mask)+im_c3(l_str.mask),repmat(clr_num_map(l_str.class),sum(sum(l_str.mask)),1)];
            end
            %[xs,ys] = ind2sub(size(im_c1),find(l_str.mask));
            %ss = [max(xs)-min(xs)+1,max(ys)-min(ys)+1];
            %im_patch = cat(3,reshape(im_c1(l_str.mask),ss),reshape(im_c2(l_str.mask),ss),reshape(im_c3(l_str.mask),ss));
            %figure, imshow(uint8(im_patch));
            %im_patch = objectAnnotation(im_patch,'Rectangle',[1,1,[size(im_patch,2)-1,size(im_patch,1)-1]],l_str.color,3);
            %figure, imshow(uint8(im_patch));
%             im_arr2{l_str.id} = objectAnnotation(double(cell2mat(im_arr2(l_str.id))),'Rectangle',round(l_str.bbox),l_str.class,7);
        end
    end
end
%labeled_data = rand_label;
%labeled_data(:,12) = labeled_data(:,12)<7;
%labeled_data(:,1:3) = labeled_data(:,1:3)/255;
rand_label=label_dataset(randperm(size(label_dataset,1)),:);
%rand_label = label_dataset;
% rr = 4;
% cc = 5;
% for ll = 1:20 
%     im = cell2mat(im_arr2(ll));
%     c = mod(ll-1,cc);
%     r = idivide(uint16(ll-1),uint16(cc));
%     im_big(r*size(im,1)+1:(r+1)*size(im,1),c*size(im,2)+1:(c+1)*size(im,2),:) = uint8(im);
% end
% imshow(im_big);
%% save the results you got
folderName = 'Results';

e_val = exist(folderName);
if ~(e_val == 7)
    mkdir(folderName);
end
cd(folderName);

c_list = ls;
b_num = sort(str2num(c_list(3:end,:)));

if isempty(b_num)
    b_num = 0;
else
    b_num = b_num(end);
end

b_num = b_num + 1;
b_str = num2str(b_num);
mkdir(b_str);
cd(b_str);

save('detectionLabels', 'rand_label');
save('Details', 'im_str');
save('Labels', 'im_labelers');
save('LabelsFolder','save_folder');
cd ../..
%%
if(version2)
    detection_test_result = ~xor(rand_label(:,12), rand_label(:,11)<7);
    %sum(rand_label(xor(rand_label(:,12),(rand_label(:,4)>0.72 & rand_label(:,10)>100)),11)~=7)
    detection_accuracy = sum(detection_test_result)/numel(detection_test_result)
else
    detection_test_result = ~xor((rand_label(:,4)>0.72 & rand_label(:,10)>50), (rand_label(:,11)<7));
    detection_accuracy = sum(detection_test_result)/numel(detection_test_result)
end

for cc = cell2mat(colors_num_class)
    cluster_data = rand_label(rand_label(:,11)==cc,:);
    detection_test_result1 = ~xor(cluster_data(:,12), cluster_data(:,11)<7);
    detection_test_result2 = ~xor((cluster_data(:,4)>0.72 & cluster_data(:,10)>50), cluster_data(:,11)<7);
    %sum(rand_label(xor(rand_label(:,12),(rand_label(:,4)>0.72 & rand_label(:,10)>100)),11)~=7)
    detection_accuracy_per_class(cc,1) = cc;
    detection_accuracy_per_class(cc,2) = sum(detection_test_result1)/numel(detection_test_result1);
    detection_accuracy_per_class(cc,3) = sum(detection_test_result2)/numel(detection_test_result2);
    %detection_accuracy_per_class(cc,4) = numel(detection_test_result1);
end
detection_accuracy_per_class
% to see whether the dataset is balanced or not...
%sum((rand_label(:,4)>0.72 & rand_label(:,10)>100) & (rand_label(:,11)<7))/numel(detection_test_result)
detection_accuracy = sum(detection_test_result)/numel(detection_test_result)