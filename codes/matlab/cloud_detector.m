function [ im ] = cloud_detector( obj, j )
%CLOUD_DETECTOR Summary of this function goes here
im = obj.imread(j);
new_im = normalize(im(:,:,1)-im(:,:,3),0,255);
figure(1000);
imshow(new_im);
pause(3);
end

