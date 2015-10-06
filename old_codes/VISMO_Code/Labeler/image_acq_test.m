im_list  = ls('FirstCameraCaptures/Images/2013_4_23/17');
im_list = cellstr(im_list(3:end,:));
im_arr = cellfun(@(x) strread(strtok(x,'.'),'%s','delimiter','_'),im_list,'UniformOutput',false);
im_to_sort = cell2mat(cellfun(@(x) str2num(x{4})+str2num(x{3})*1000+str2num(x{2})*100000+str2num(x{1})*10000000,im_arr,'UniformOutput',false));
[~,sorted_indices]=sort(im_to_sort);

im_list = im_list(sorted_indices);
cd 'FirstCameraCaptures/Images/2013_4_23/17';
for i = 1:length(im_list)/3
    cim_list{i} = uint8(double(imread(im_list{((i-1)*3)+1}))/3 ...
        + double(imread(im_list{((i-1)*3)+2}))/3 ...
        + double(imread(im_list{((i-1)*3)+3}))/3);
end

cd ../../../..
LabelerGUI(cim_list)
