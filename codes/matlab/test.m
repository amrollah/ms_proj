% close all;
% clear all;
% pm=6000;
% addpath(genpath('C:\Users\chamsei\Documents\GitHub\ms_proj\old_codes\vismolib'))
% vmlconfig_cavriglia;
% conf = evalin('base','VMLCONF');
% obj = vmlSeq('2015_08_10',[6 20]);
% 
% ts = timeseries(obj.Irr(:,2:end), obj.Irr(:,1));
% ts2 = resample(ts, obj.P(:,1)); 
% figure;
% plot(ts2.data(1:end-pm,1), obj.P(1:end-pm,2),'b.');
% hold on;
% plot(ts2.data(end-pm:end,1), obj.P(end-pm:end,2),'r.');
% 
% 
% degree = 130;
% theta = degtorad(degree); % rotaion angle counterclockwise about the origin
% R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% theta_inv = degtorad(-degree); % rotaion angle counterclockwise about the origin
% R_inv = [cos(theta_inv) -sin(theta_inv); sin(theta_inv) cos(theta_inv)];
% 
% translationVector = [0;0];
% 
% figure;
% plot(obj.Irr(obj.clear_times,1),obj.Irr(obj.clear_times,2), '.g');
% hold on;
% plot(obj.Irr(obj.cloudy_times,1),obj.Irr(obj.cloudy_times,2), '.r');
% 
% figure;
% load('all_sun_pos__2015_09_21.mat');
% plot(all_pos(:,1),all_pos(:,2), '-r');
% hold on;
% plot(all_pos(:,3),all_pos(:,4), '-g');
% 
% % create a matrix which will be used later in calculations
% center_pos = fix((all_pos(1,1:2)+all_pos(end,1:2))/2); %all_pos(fix(size(all_pos,1)/2),1);
% [center_pos, ~] = sunpos_midday(obj);
% %y_center_pos = all_pos(fix(size(all_pos,1)/2),2);
% center_pos = repmat(center_pos, 1, size(all_pos,1));
% 
% % create a matrix which will be used later in calculations
% center_impos = fix((all_pos(4,3:4)+all_pos(end-4,3:4))/2);
% center_impos = repmat(center_impos', 1, size(all_pos,1));
% 
% 
% % shift points in the plane so that the center of rotation is at the origin
% rot_pos = R*(all_pos(:,1:2)' - center_pos) + center_pos;
% rot_impos = R*(all_pos(:,3:4)' - center_pos) + center_pos;
% 
% sh = .030;
% [min_y, idx] = min(rot_pos(2,:));
% min_x = rot_pos(1,idx);
% rot_pos(2,:) = (rot_pos(2,:)-min_y).*(1+(sh^2)) + min_y;
% rot_pos(1,:) = (rot_pos(1,:)-min_x).*(1-sh) + min_x;
% %rot_pos(:,:) = (rot_pos(1:2,:)-[min_x;min_y]).*[1-sh;1+(sh^2)] + [min_x;min_y];
% 
% inv_center_pos_pt = fix((rot_pos(1:2,1)+rot_pos(1:2,end))/2);
% inv_center_pos = repmat(inv_center_pos_pt, 1, size(rot_pos,2));
% orig_pos = R_inv*(rot_pos - center_pos) + center_pos;
% 
% figure;
% plot(rot_pos(1,:),rot_pos(2,:), '-r');
% hold on;
% plot(rot_impos(1,:),rot_impos(2,:), '-b');
% 
% figure;
% 
% plot(all_pos(:,1),all_pos(:,2), '-r');
% hold on;
% plot(all_pos(:,3),all_pos(:,4), '-g');
% hold on;
% plot(orig_pos(1,:),orig_pos(2,:), '-b');
% 
% for j=1:length(cloudy)
%     disp(datevec(obj.ti(cloudy(j))));
% end


% figHandles = get(0,'Children'); % gets the handles to all open figure
% rez=800; %resolution (dpi) of final graphic
% resolution=get(0,'ScreenPixelsPerInch'); %dont need to change anything here
% 
% for f = figHandles'
% %     axesHandle = get(f,'Children'); % get the axes for each figure
% %     titleHandle = get(axesHandle,'Title'); % gets the title for the first (or only) axes)
% %     text = get(titleHandle,'String'); % gets the text in the title
% %     text = strrep(text, '/', '_');
%     text = 'sun_samples';
%     figpos=getpixelposition(f); %dont need to change anything here    
%     set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); %dont need to change anything here
%     print(f,strcat('C:\data\cav\',text),'-djpeg',['-r',num2str(rez)],'-opengl') %save file 
% 
%     %saveas(f, strcat('C:\data\cav\log_data\power_plots\',text), 'jpg') % save the figure with the title as the file name 
% end
% coef=90000/70450;
% for i=1:length(data)
%     disp(data{i}.clouds);
% %     data{i}.clouds=data{i}.clouds*coef;
% end

% % load('D:\abb_data\cc_data.mat','cc');
% for rr=1:length(rg)
%     ccp=cc(rg(rr)+1:rg(rr+1));
%     save(['D:\abb_data\cc_data',num2str(rr), '.mat'], 'ccp');
% end
% % load('D:\abb_data\cc_data.mat','cc');
% cc_all = cellfun(@(d) sum(d.cc_fact), data);
% clouds = cellfun(@(d) d.clouds, data);
% time = cellfun(@(d) d.time, data);
% TV  = datevec(time);
% mm=mean(TV(inx,4:5));
% inx = find(isnan(cc_all));
% img_save_path = 'E:\ABB\cav\diffuse\';
% img_save_path2 = 'E:\ABB\cav\diffuse_cc_0\';
% for pp=1:length(inx)
%     disp(pp);
%     d=data{inx(pp)};
%     nan_idx=find(isnan(d.cc_fact));
%     d.cc_fact(nan_idx)=100;
%     data{inx(pp)}=d;
% %     disp(tm(4:5));
% %     filename=[d.day '__' num2str(d.j) '.jpeg'];
% %     copyfile([img_save_path filename],strcat(img_save_path2,filename));  
% end
% close(1);
figure(1);
data1 = randn(1,1e5); %// example data
data2 = randn(1,1e5) + .5*data1; %// example data correlated to above
values = hist3([data1(:) data2(:)],[51 51]);
imagesc(values)
newmap = jet;
newmap(1,:) = [1 1 1];
colormap(newmap);    
colorbar
axis equal
axis xy