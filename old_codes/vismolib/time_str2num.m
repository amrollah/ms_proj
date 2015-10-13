function t = time_str2num(s)
y=textscan(s,'%s','delimiter',':.');
y = y{1};
if length(y)~=2, t=NaN; else t = [1 1/60]*cellfun(@str2double,y); end
