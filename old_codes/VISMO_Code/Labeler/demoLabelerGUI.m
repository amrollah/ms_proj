
%%
%im_dir = 'imagesets/2008/';
%im_dir = 'C:/Users/chbuzey/Documents/Visual Studio 2010/Projects/ImageAcquisitionTool/ImageAcquisitionTool/bin/Debug/Images/2013/05/14/';
im_dir = 'C:\Users\chbuzey\Documents\MATLAB\Archive\Validation Sequences\Cloud Detection\';  % for the validation of detection
%im_dir = 'C:\Users\chbuzey\Documents\MATLAB\Archive\Validation Sequences\Cloud Tracking\Sequence2\';  % for the validation of detection

rect = 0;
e_val = exist('im_arr');

tic
if(rect)
    e_val1 = exist('lut');
    if(e_val1~=1)
        lut=getLUT(double(imread('sky_image0.jpg')),1);
    end
    if(e_val ~= 1)
        im_list  = ls(im_dir);
        im_list = cellstr(im_list(3:end,:));
        im_arr = cellfun(@(x) sph2rect(double(imread([im_dir x])),lut),im_list(:),'UniformOutput',false);
    end
else
    if(e_val ~= 1)
        im_list  = ls(im_dir);
        im_list = cellstr(im_list(3:end,:));
        im_arr = cellfun(@(x) imread([im_dir x]),im_list,'UniformOutput',false);
    end
end
toc

details_struct = struct('im_dir',{im_dir},'im_list',{im_list});
labeler_struct = struct('im_arr',{im_arr},'im_str',{details_struct});

%LabelerGUI(labeler_struct,im_labelers);
LabelerGUI(LabelerPrep('details',im_str),im_labelers);

%%
% any labeled pixel has the following set of attributes: 
% red, green, blue, red/blue, red/green, green/blue, red-blue, red-green,
% green-blue, red+blue+green, label
colors_num_class = {1,2,3,4,5,6,7,8,9};
colors_num_keys = {'m', 'g', 'y', 'b', 'p', 'r', 'k', 'w', 'c'};
clr_num_map = containers.Map(colors_num_keys(:),colors_num_class(:));
label_dataset = [];
for i = length(im_labelers):-1:1
    for j = 1:length(im_labelers(i).im_labels)
        
        l_str = im_labelers(i).im_labels(j);
        im_c = cell2mat(im_arr(l_str.id));
        im_c1 = double(im_c(:,:,1));
        im_c2 = double(im_c(:,:,2));
        im_c3 = double(im_c(:,:,3));
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
        [xs,ys] = ind2sub(size(im_c1),find(l_str.mask));
        ss = [max(xs)-min(xs)+1,max(ys)-min(ys)+1];
        im_patch = cat(3,reshape(im_c1(l_str.mask),ss),reshape(im_c2(l_str.mask),ss),reshape(im_c3(l_str.mask),ss));
        %figure, imshow(uint8(im_patch));
        im_patch = objectAnnotation(im_patch,'Rectangle',[1,1,[size(im_patch,2)-1,size(im_patch,1)-1]],l_str.color,3);
        %figure, imshow(uint8(im_patch));
    end
end

rand_label=label_dataset(randperm(size(label_dataset,1)),:);

%csvwrite('labeled_pixels_cloud_detection.csv',rand_label);