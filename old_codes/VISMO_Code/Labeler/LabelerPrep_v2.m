function [ labeler_struct_v2 ] = LabelerPrep_v2( varargin )

vargs = varargin;
nargs = length(vargs);
names = vargs(1:2:nargs);
values = vargs(2:2:nargs);

validnames = {'dir', 'details'};    
for name = names
   validatestring(name{:}, validnames);
end

match_s = strmatch('dir',names);
if(~isempty(match_s))
    im_dir = values{match_s};
    im_list = ls(im_dir);
    im_list = cellstr(im_list(3:end,:));
    details_struct = struct('im_dir',{im_dir},'im_list',{im_list});
end

match_d = strmatch('details', names);
if(~isempty(match_d))
    details_struct = values{match_d};
    im_dir = details_struct.im_dir;
    im_list = details_struct.im_list;
end

match_d = strmatch('rec', names);
if(~isempty(match_d))
    recordings_folder = values{match_d};
    details_struct.rec_folder = recodings_folder;
end

if(~isempty(match_d) && ~isempty(match_s))
    warning('Should not specify both the directory and the struct at the same time!!');
end
%im_arr = cellfun(@(x) imread([im_dir x]),im_list,'UniformOutput',false);
%labeler_struct = struct('im_arr',{im_arr},'im_str',{details_struct});
labeler_struct_v2 = details_struct;
end

