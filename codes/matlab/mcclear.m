close all; clear all; clc;

%url = 'http://www.soda-pro.com/portlets-common/cgi-bin/proxy.py?url=http://toolbox.webservice-energy.org/service/wps';
filename = 'C:\data\mcclear\irradiation-clear_sky_cav_2014.csv';
fid = fopen(filename);
out = textscan(fid,'%s%f%f%f%f%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s','delimiter',';');
fclose(fid);

 rmse=sqrt(sum((data(:)-estimate(:)).^2)/numel(data));
 