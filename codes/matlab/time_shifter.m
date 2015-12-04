function [ hour,min,sec ] = time_shifter(obj,tm)
  [year,month,day,hour,min,sec]=datevec(tm);
  % for Daylight saving
  if ~isdst(datetime(year,month,day,hour,min,sec,'TimeZone',obj.conf.timezone))
    hour = hour + 1;
  end
end

