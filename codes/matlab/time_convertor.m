function [ str ] = time_convertor( t )
%TIME_CONVERTOR converts time in milliseconds to readable format

str = strcat(num2str(t(4)),':',num2str(num2str(t(5)))); %,':',num2str(num2str(floor(t(6)))

end

