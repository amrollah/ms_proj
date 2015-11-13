function [ t ] = time_struct(obj,t_in)
%CONSTRUCT_TIME converts a hour:min time format to full datetime format
parts = strsplit(char(t_in),':');
t = [];
[t.year,t.month,t.day]=datevec(obj.ti(1));
t.hour = str2double(parts(1));
t.minute = str2double(parts(2));
t.second = 0;
t.UTCOffset = obj.calib.model3D.UTC;

end

