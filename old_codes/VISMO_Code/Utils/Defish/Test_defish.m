%% Test Script for Fish Eye Defishing Code with 9-Point interpolation.
%
% B. Zeydan, 05. Mar. 2013

[sky_im] = cv.imread('lowest.jpg');
%[sky_im] = cv.imread('ChessBoardImage.jpg');
%[sky_im] = cv.imread('calib1.jpg');
sky_im = double(sky_im);
%sky_im = cv.resize(sky_im,[500,500]);
figure, imshow(uint8(sky_im));
%%
%def_im = sph2rect(sky_im);
%figure, imshow(def_im);
%%

if(exist('lut')~=1)
    lut = getLUT(sky_im,1,ocam_model);
end
tic
def_im = sph2rect(sky_im,lut);
toc
%cv.imwrite('fixed.jpg',def_im);
figure, imshow(uint8(def_im));