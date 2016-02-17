function t = time_str2num(s)
if isempty(s), t=NaN; return; end
y=textscan(s,'%s','delimiter',':.');
y = cellfun(@str2double,y{1});
f = [1 1/60 1/3600]/24;
if length(y)>3, t=NaN; else t = f(1:length(y))*y; end
