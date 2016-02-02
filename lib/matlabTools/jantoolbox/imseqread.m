function imseq = imseqread(filter)
folder = [fileparts(filter) filesep];
files = dir(filter);
if isempty(files), error(['no image matching filter ' filter]); end
imseq = cell(1,length(files));
for i=1:length(files)
  imseq{i}=imread([folder files(i).name]);
end
