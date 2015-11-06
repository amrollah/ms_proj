%close all;
%% read McClear data file (estimated data)
mcclear_path = 'C:\data\mcclear\';
filename_mcclear = strcat(mcclear_path, 'irradiation-clear_sky_cav_2014_v2.csv');
fid = fopen(filename_mcclear);
mc_data = textscan(fid,'%s%s%f%f%f%f%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s','delimiter',';');
fclose(fid);
i=10;


for i=10:24
    try
        figure(i);
        plot(mc_data{1,i}(721:1440:size(mc_data{1,1},1)), '.');
    catch
        close(i);
    end
end
ColorSet = varycolor(365);

figure(20);
hold all;
for i=0:365
    plot(mc_data{1,4}((1440*i)+1:1440*(i+1))*60, '-', 'Color', ColorSet(i+1,:));
    pause(1);
end

t = abs(mc_data{1,4}(721:1440:(size(mc_data{1,1},1)-1440)) - mc_data{1,4}(721+1440:1440:size(mc_data{1,1},1)));

filename_mcclear = strcat(mcclear_path, 'irradiation-clear_sky_cav_2013_v2.csv');
fid = fopen(filename_mcclear);
mc_data_2013 = textscan(fid,'%s%s%f%f%f%f%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s','delimiter',';');
fclose(fid);

filename_mcclear = strcat(mcclear_path, 'irradiation-clear_sky_cav_2014_v2.csv');
fid = fopen(filename_mcclear);
mc_data_2014 = textscan(fid,'%s%s%f%f%f%f%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s','delimiter',';');
fclose(fid);

filename_mcclear = strcat(mcclear_path, 'irradiation-clear_sky_cav_2015_till_september_v2.csv');
fid = fopen(filename_mcclear);
mc_data_2015 = textscan(fid,'%s%s%f%f%f%f%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s','delimiter',';');
fclose(fid);

day_idx = mc_data_2013{1,5}>1;
diffuse = mc_data_2013{1,6}(day_idx);
direct = mc_data_2013{1,5}(day_idx);
r_2013 = (diffuse./direct)*100;
diff_2013 = abs(r_2013(1:end-1) - r_2013(2:end));


day_idx2 = mc_data_2014{1,5}>1;
diffuse2 = mc_data_2014{1,6}(day_idx2);
direct2 = mc_data_2014{1,5}(day_idx2);
r_2014 = (diffuse2./direct2)*100;
diff_2014 = abs(r_2014(1:end-1) - r_2014(2:end));