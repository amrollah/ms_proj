close all;
clear all;
clc;
training = false;
img_save_path='';
prj_path='';
proj_path;

addpath(prj_path);
s = vml('cavriglia','2015_08_03',[],true);
conf = [];
conf=local_conf(conf);
load('calc\data_clean.mat', 'data');

if training   
    imgs = dir(img_save_path); % get all content of Cavriglia image folder
    imgs = imgs(~vertcat(imgs.isdir)); % filter only folders
    sat_vals = {};
    all_patches = [];

    for i=1:length(data)
        d=data{i};
        if d.sun_flag==4
            continue;
        end
        img_file = [img_save_path d.day '__' num2str(d.j) '.jpeg'];
        I=imread(img_file);
        sun_pos = s.sunpos_im(d.time);
        cx=floor(sun_pos(1)); cy=floor(sun_pos(2)); R=100;
        sun_patch = I(max(cx-R,1):min(cx+R,size(I,1)),max(cy-R,1):min(cy+R,size(I,2)),:); 
        h=figure(1);
    %     maxfig(h,1);
        [x, y, BW, xi, yi] = roipoly(sun_patch);
        rI = reshape(sun_patch,[],3);
        patch = rI(BW,:);
        all_patches = [all_patches; patch];
    %     figure(3); imshow(BW);
        pause();
    end
    all_patches=double(all_patches);
    save('calc\all_patches.mat','all_patches');
    figure(5);
    subplot(2,2,1);
    rh=histogram(all_patches(:,1),10);
    title('red');
    subplot(2,2,2);
    gh=histogram(all_patches(:,2),10);
    title('green');
    subplot(2,2,3);
    bh=histogram(all_patches(:,3),10);
    title('blue');

else
    load('calc\all_patches.mat','all_patches');
    [rh,rc]=hist(all_patches(:,1),10);
    [gh,gc]=hist(all_patches(:,2),10);
    [bh,bc]=hist(all_patches(:,3),10);
    p_w = 4;
    R = 100;
    thres=200;
    %evaluation
    for i=1:length(data)
        d=data{i};
        if d.sun_flag==4
            continue;
        end
        img_file = [img_save_path d.day '__' num2str(d.j) '.jpeg'];
        I=imread(img_file);
        sun_pos = s.sunpos_im(d.time);
        cx=floor(sun_pos(1)); cy=floor(sun_pos(2));
        xgrid = max(cx-R,1):p_w:min(cx+R,size(I,1));
        ygrid = max(cy-R,1):p_w:min(cy+R,size(I,2));
        [X,Y] = meshgrid(xgrid,ygrid);
%         h=figure(10);
%         imshow(I);hold on;
        figure(20);
        for p=1:numel(X)
            patch = imcrop(I,[X(p),Y(p),p_w-1,p_w-1]);
            subplot(ceil(sqrt(numel(X))),ceil(sqrt(numel(X))),p);
            imshow(patch);
            patch=reshape(patch,[],3);
            prh=hist(patch(:,1),rc);
            pgh=hist(patch(:,2),gc);
            pbh=hist(patch(:,3),bc);
            distance = mean([pdist2(prh,rh),pdist2(pgh,gh),pdist2(pbh,bh)]);
            disp(distance);
            if distance<thres
%                 rectangle('Position',[X(p),Y(p),p_w-1,p_w-1]);
            end
        end
        pause();
    end
end