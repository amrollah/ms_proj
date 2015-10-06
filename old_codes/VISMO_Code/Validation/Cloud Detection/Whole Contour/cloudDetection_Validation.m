function [h, class] = cloudDetection_Validation(varargin)
% Visualizes the polygons in an image.
%
% LMplot(annotation, img)
% or
% LMplot(database, ndx, HOMEIMAGES)
%
% Example:
%   [annotation, img] = LMread(filename, HOMEIMAGES)
%   LMplot(annotation, img)
%
%   thumbnail = size of the image
%   The plot uses only 7 colors, therefore, some times, different object.
%   Classes might have the same color assigned.
%
% If an object has the field 'confidence', the thickness of the bounding
% box will be equal to the confidence.
%
% If the object has the filed 'detection', then when the state is 'false'
% it will show the outline in green, and when the state is 'correct' it
% will show the outline in red.

switch length(varargin)
   case 1
       annotation = varargin{1};
       img = -1;
   case 2
       annotation = varargin{1};
       img = varargin{2};
   case 3
       D = varargin{1};
       ndx = varargin{2};
       
       if ndx>length(D)
           return
       end
       
       HOMEIMAGES = varargin{3};

       annotation = D(ndx).annotation;
       HOMEIMAGES
       %annotation.folder = 'users\burakzeydan\cloud_images';
       img = imread(fullfile(HOMEIMAGES, 'users\burakzeydan\cloud_images', annotation.filename));
       %img = LMimread(D, ndx, HOMEIMAGES); % Load image
   otherwise
       error('Too many input arguments.')
end

% Define colors
colors = 'rgbcmyw';
colors = hsv(15);

% Draw each object (only non deleted ones)
h = []; class = [];

if isfield(annotation, 'object')
   Nobjects = length(annotation.object); n=0;
   for i = 1:Nobjects
       n = n+1;
       %class{n} = annotation.object(i).name; % get object name
       %col = colors(mod(sum(double(class{n})),15)+1, :);
       [X,Y] = getLMpolygon(annotation.object(i).polygon);
       X = double(X(:)); Y = double(Y(:));
       contours{1,i} = num2cell([X,Y],2)';
   end
   im_labeled = zeros(size(img));
   im_labeled = cv.drawContours(im_labeled,contours,'Thickness',-1);
   figure, imshow(im_labeled);
   
   sunPattern = imread('sunPattern.jpg');
   load 'image_mask2.mat';
   img(~image_mask) = 0;  % we eliminate the parts that are angularly behind the camera observation point
   
   res = cv.matchTemplate(img,sunPattern);
   [x,y]=ind2sub(size(res),find(res==min(min(res))));
   img = cv.circle(img,[y+3,x+1],130,'Color',[0,0,0],'Thickness',-1);
   
   r_b = double(img(:,:,1))./double(img(:,:,3)); % R/B thresholding applied
   rgb = double(img(:,:,1))+double(img(:,:,2))+double(img(:,:,3)); % to distinguish other objects like buildings from clouds, a test on brightness
   mask = r_b > 0.72 & rgb > 100;             %% brightness auto thing is needed because of this ..
   figure, imshow(mask);
   
end

