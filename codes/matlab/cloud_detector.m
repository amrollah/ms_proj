function [ im, contours, bboxes ] = cloud_detector(obj, j, sunPosition)
%CLOUD_DETECTOR: detects cloud patches in an image
im = obj.imread(j);
load([obj.conf.datafolder obj.conf.calibration{1}]);
load([obj.conf.datafolder 'calib.mat']);

%obj.conf.ocam_model = ocam_model;
%new_im = normalize(im(:,:,1)-im(:,:,3),0,255);
maskMethod = 'body';
sunPatch = getSunPatch(obj,im,sunPosition);
im = eliminateSun(obj,im,sunPosition);
mask = getBinaryImage(im, maskMethod, sunPatch);
imshow(mask);
pause(5);
[contours, contoursProj, areas, areasProj, centroids, centroidsProj, bboxes, bboxesProj] = detectObjects(obj,mask,ocam_model,R);
%mask = objectAnnotation(mask,'Rectangle',bboxes,'g',3);
%imshow(mask);
im = cv.drawContours(im,contours,'Color',[255 255 0],'Thickness',3);
%figure(1000);
%imshow(im);
%pause(1);
end
%%
function sunPatch = getSunPatch(obj, im, sunPosition)
    [X,Y] = meshgrid(1:size(im,1),1:size(im,2));
    sunPatch = sqrt((X-sunPosition(2)).^2+(Y-sunPosition(1)).^2)<=obj.conf.circumsolarSize; % sun patch for better binarization is calculated according to the sun position on the image
end
%%
function [frame] = eliminateSun(obj,frame,sunPosition)
    x = sunPosition(2);
    y = sunPosition(1);        
    frame = cv.circle(frame,[x,y],obj.conf.circumsolarSize,'Color',[0,0,0],'Thickness',-1); % the sun is eliminated using a patch
end
%%    
function [contours, contoursProj, areas, areasProj, centroids, centroidsProj, bboxes, bboxesProj] = detectObjects(obj,mask,ocam_model,R)
    %contours = cv.findContours(mask,'Mode','List'); % contour extraction using openCV method
    contours = cv.findContours(mask,'Mode','External'); % contour extraction using openCV method

    % areas of objects are approximately computed here using several methods
    %areas = cell2mat(cellfun(@(x) size(cell2mat(x(:)),1),contours,'UniformOutput',false))';
    %areas = cell2mat(cellfun(@(x) prod(max(cell2mat(x(:)),[],1)-min(cell2mat(x(:)),[],1)) ,contours,'UniformOutput',false))';
    areas = cell2mat(cellfun(@(x) cv.contourArea(x(:)) ,contours,'UniformOutput',false))';

    % elminate the objects with contour smaller than 350.
    consider = (areas > 600); % if the area of an object is less than 300 pixels, we discard that

    % take only the objects that we want to consider
    contours = contours(consider);
    areas = areas(consider);

    % contours are projected onto the hypothetical plane at an height
    % of cbh in the world coordinate system then the area of the
    % projected contours are computed.
    contoursProj = cellfun(@(x) num2cell(getImage2Real(cell2mat(x(:))+1,ocam_model,R,obj.conf.cbh),2)',contours,'UniformOutput',false);
    areasProj = cell2mat(cellfun(@(x) cv.contourArea(x(:)) ,contoursProj,'UniformOutput',false))';

    % different methods for determining the point that represent the
    % cloud objets first method uses the mean of the contours the
    % second one uses the mid point of the bounding rectangle
    %centroids = cellfun(@(x) mean(cell2mat(x(:)),1)+1,contours,'UniformOutput',false);
    %centroidsProj = cellfun(@(x) mean(cell2mat(x(:)),1),contoursProj,'UniformOutput',false);
    centroids = cellfun(@(x) (min(cell2mat(x(:))+1,[],1)+max(cell2mat(x(:))+1,[],1))/2,contours,'UniformOutput',false);
    centroidsProj = cellfun(@(x) (min(cell2mat(x(:)),[],1)+max(cell2mat(x(:)),[],1))/2,contoursProj,'UniformOutput',false);

    % bounding boxes on the image and the projected plane are
    % deteremined
    bboxes = cellfun(@(x) [min(cell2mat(x(:))+1,[],1),max(cell2mat(x(:))+1,[],1)-min(cell2mat(x(:))+1,[],1)] ,contours,'UniformOutput',false);
    bboxesProj = cellfun(@(x) [min(cell2mat(x(:)),[],1),max(cell2mat(x(:)),[],1)-min(cell2mat(x(:)),[],1)] ,contoursProj,'UniformOutput',false);

    % centroids and boxes on both planes are put into 2D arrays for
    % easy access later
    centroids = vertcat(centroids{:});
    centroidsProj = vertcat(centroidsProj{:});
    bboxes = vertcat(bboxes{:});
    bboxesProj = vertcat(bboxesProj{:});
end
%%
function mask = getBinaryImage(frame, type, varargin)
    % two types are available: 'body' and 'edge'.
    r_b = double(frame(:,:,1))./double(frame(:,:,3)); % R/B thresholding applied
    rgb = double(frame(:,:,1))+double(frame(:,:,2))+double(frame(:,:,3)); % to distinguish other objects like buildings from clouds, a test on brightness
    mask = r_b > 0.92 & rgb > 250;             %% brightness auto thing is needed because of this ..

    % this part is for the sun glare elimination. 
    if(~isempty(varargin))
        sunP = varargin{1}; % the area around the sun is defined by the sun patch
        g_b = double(frame(:,:,2))./double(frame(:,:,3)); % G/B thresholding is applied for the parts around the sun
        mask2 = g_b > 0.96;  % apply another threshold on the G/B feature.
        %mask2 = (r_b - g_b)>-0.01;
        mask(sunP) = mask2(sunP);
    end

    % if the cloud binarization is done using edges, this extra part is
    % run
    if(strcmp(type,'edge')) 
        H = [1,1,1;1,-8,1;1,1,1];
        mask = imfilter(mask,H,'replicate');
    end
    mask = logical(cv.medianBlur(mask, 'KSize', 13)); % the binary image is filtered to eliminate the possible noise
end
%%
function pointReal = getImage2Real(m,calibModel,R,cbh)
    % here we project a set of points given in m onto a plane that is cbh meters away from the camera the used convention in opencv is [column,row] kind of
    % represantation for pixel points, here however we change that to [row,column]
    R_earth = 6378100; 
    pointTemp = cam2world([m(:,2),m(:,1)]',calibModel); % we get the corresponding direction on a unit-sphere at first.
    pointReal = zeros(size(pointTemp,2),2); % we crate the array that will contain the projected positions
    if(size(R,1)==3)
        pointTemp = R*pointTemp;
        alphas = (pi/2) - atan2(-pointTemp(3,:),sqrt(pointTemp(1,:).^2 + pointTemp(2,:).^2));
        heights = (-2*R_earth*cos(alphas) + sqrt((2*R_earth*cos(alphas)).^2 + 4*(cbh^2+2*cbh*R_earth)))/2;
        pointReal = (pointTemp(1:2,:).*repmat(heights,2,1))';
    else
        pointReal = ((R*pointTemp(1:2,:)).*repmat((-cbh./pointTemp(3,:)),2,1))'; % here we map the point on the unit-sphere onto the plane
    end
end
