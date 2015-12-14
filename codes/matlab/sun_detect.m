clear all;
close all;
clc;

load('sun_patches.mat');
max_r=0;
max_g=0;
max_b=0;
for i=1:size(sun_patches,2)
    figure(1);
    subplot(14,15,i);
    rgb = sun_patches{1,i};
    imshow(rgb);
    sun=imcrop(rgb,sun_patches{2,i});
    
    for d=1:3
    figure(d+1);
    subplot(14,15,i);    
    [f,x]=hist(double(reshape(sun(:,:,d),1,[])),10);
    bar(x,f/sum(f));
    axis off;
    xlim([0,260]);
    ylim([0,1]);
    end
%     max_r = max(max_r,max(max(sun(:,:,1))));
%     max_g = max(max_g,max(max(sun(:,:,2))));
%     max_b = max(max_b,max(max(sun(:,:,3))));
    pause(0.1);
end