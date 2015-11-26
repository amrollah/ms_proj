close all;
clear all;
degree = 130;
theta = degtorad(degree); % rotaion angle counterclockwise about the origin
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
theta_inv = degtorad(-degree); % rotaion angle counterclockwise about the origin
R_inv = [cos(theta_inv) -sin(theta_inv); sin(theta_inv) cos(theta_inv)];

translationVector = [0;0];

vmlconfig_cavriglia;
conf = evalin('base','VMLCONF');
obj = vmlSeq('2015_08_21',[6 20]);

figure;
plot(obj.Irr(obj.clear_times,1),obj.Irr(obj.clear_times,2), '.g');
hold on;
plot(obj.Irr(obj.cloudy_times,1),obj.Irr(obj.cloudy_times,2), '.r');

figure;
load('all_sun_pos__2015_09_21.mat');
plot(all_pos(:,1),all_pos(:,2), '-r');
hold on;
plot(all_pos(:,3),all_pos(:,4), '-g');

% create a matrix which will be used later in calculations
center_pos = fix((all_pos(1,1:2)+all_pos(end,1:2))/2); %all_pos(fix(size(all_pos,1)/2),1);
[center_pos, ~] = sunpos_midday(obj);
%y_center_pos = all_pos(fix(size(all_pos,1)/2),2);
center_pos = repmat(center_pos, 1, size(all_pos,1));

% create a matrix which will be used later in calculations
center_impos = fix((all_pos(4,3:4)+all_pos(end-4,3:4))/2);
center_impos = repmat(center_impos', 1, size(all_pos,1));


% shift points in the plane so that the center of rotation is at the origin
rot_pos = R*(all_pos(:,1:2)' - center_pos) + center_pos;
rot_impos = R*(all_pos(:,3:4)' - center_pos) + center_pos;

sh = .030;
[min_y, idx] = min(rot_pos(2,:));
min_x = rot_pos(1,idx);
rot_pos(2,:) = (rot_pos(2,:)-min_y).*(1+(sh^2)) + min_y;
rot_pos(1,:) = (rot_pos(1,:)-min_x).*(1-sh) + min_x;
%rot_pos(:,:) = (rot_pos(1:2,:)-[min_x;min_y]).*[1-sh;1+(sh^2)] + [min_x;min_y];

inv_center_pos_pt = fix((rot_pos(1:2,1)+rot_pos(1:2,end))/2);
inv_center_pos = repmat(inv_center_pos_pt, 1, size(rot_pos,2));
orig_pos = R_inv*(rot_pos - center_pos) + center_pos;

figure;
plot(rot_pos(1,:),rot_pos(2,:), '-r');
hold on;
plot(rot_impos(1,:),rot_impos(2,:), '-b');

figure;

plot(all_pos(:,1),all_pos(:,2), '-r');
hold on;
plot(all_pos(:,3),all_pos(:,4), '-g');
hold on;
plot(orig_pos(1,:),orig_pos(2,:), '-b');

for j=1:length(cloudy)
    disp(datevec(obj.ti(cloudy(j))));
end
