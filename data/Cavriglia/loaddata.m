folder1 = textscan(folder,'%f','delimiter','_');
folder1 = strjoin(cellfun(@num2str,num2cell(reshape(folder1{1}(3:-1:1),1,3)),'UniformOutput',0),'_');
tm1 = datenum(folder,'yyyy_mm_dd');

idstr = [];
x1 = dir([conf.basefolder folder '\*.7z']);
for i=1:length(x1)
  j = strfind(x1(i).name,folder1);
  if ~isempty(j), idstr = x1(i).name(12:j-2); end
  cmd = ['"C:\Program Files\7-Zip\7z.exe" x -y -o' conf.basefolder folder ' ' conf.basefolder folder '\' x1(i).name];
  disp(cmd);
  system(cmd);
end

datalogfile = [conf.basefolder folder '\' folder '_data.log'];
if exist(datalogfile,'file')
  disp(['READING FROM DATA LOGGER FILE ' datalogfile]);
  %[x1,tx1] = swallowcsv(datalogfile,' ');
  fid = fopen(datalogfile);
  x1 = textscan(fid, '%s%s%f%f%f%f%f%f','delimiter',' ');
  fclose(fid);
  tx1_times=cellfun(@(s)textscan(s,'%f','delimiter',':'),x1{1,2});
  tm = tm1+([3600 60 1]*cat(2,tx1_times{:})/86400)';
  data_Irr = [tm x1{1,5} x1{1,6}];
  data_Temp = [tm x1{1,7} x1{1,8}];
else 
  disp(['MISSING DATA LOGGER FILE ' datalogfile]);
  data_Irr = [];
  data_Temp = [];
end

powerfile = [conf.basefolder folder '\' idstr '-' folder1 '.txt'];
if ~isempty(idstr) && exist(powerfile,'file')
  disp(['READING FROM POWER LOG FILE ' powerfile]);
  fid = fopen(powerfile);
  x1 = textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',' ');
  fclose(fid);
  x1 = cell2mat(x1);
  %x1 = swallowcsv(powerfile,char(9));
  tm = tm1+x1(:,1)/86400;
  pp = x1(:,[2 12 22 32]);
  psum = sum(pp,2);
  jj=sum(pp==0,2);
  jj = ([true;jj(1:end-1)<jj(2:end)] | [true(2,1);jj(1:end-2)<jj(3:end)]) & ...
       ([jj(2:end)<jj(1:end-1);true] | [jj(3:end)<jj(1:end-2);true(2,1)]);
  psum(jj) = max(0,interp1(tm(~jj),psum(~jj),tm(jj),'linear','extrap'));
  data_Power = [tm psum pp all(x1(:,[5 15 25 35])==4,2)];
else
  disp(['MISSING POWER LOG FILE ' powerfile]);
  data_Power = [];  
end
