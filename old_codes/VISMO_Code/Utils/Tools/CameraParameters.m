%% Camera filters and conversion factor implementation.
% xsensor   : size of the long edge of the camera sensor in m.
% ysensor   : size of the short edge of the camera sensor in m.
% rsensor   : size of the diagonal of the camera sensor in m.
% refsensor : size of the radius of the image on the sensor in m.
% rpixel    : size of a pixel in m.
% cbh       : estimated cloud base height for a frame in m.
% fov       : effective fov of the lens for the sensor in degrees.
% Camera parameter and cbh dependent conversion factor calculations. Can
% become a seperate file...
% B. Zeydan, 27. Feb. 2013
%
%
% The function is seperated from the main code.
% Burak Zeydan, 22. Aug. 2014
    function camParams = CameraParameters()
        
        xsensor = 0.0034;
        xefsensor = 0.0034;
        
        ysensor = 0.0019;
        yefsensor = 0.0019;
        
        rsensor = sqrt((xsensor/2)^2 + (ysensor/2)^2); %this depends on the camera sensor..
        refsensor = sqrt((xefsensor/2)^2 + (yefsensor/2)^2); %this depends on we utilize the camera sensor..
        
        frameRes = size(frame);
        
        rpixel = xsensor/frameRes(2);
        
        fov = 90 * 0.0174532925;
        
        cbh_l = 500;
        
        timeConversionFactor = 1;
        lengthConversionFactor = cbh_l / sqrt(refsensor^2  + ((xefsensor/2)/tan(fov/2))^2);
        conversionFactor = 1;%rpixel * lengthConversionFactor / timeConversionFactor;
        
        camParams = struct(...
            'xsensor', xsensor, ...
            'ysensor', ysensor, ...
            'xefsensor', xefsensor, ...
            'yefsensor', yefsensor, ...
            'rsensor', rsensor, ...
            'refsensor', refsensor, ...
            'rpixel', rpixel, ...
            'lengthConversionFactor', lengthConversionFactor, ...
            'timeConversionFactor', timeConversionFactor, ...
            'conversionFactor', conversionFactor, ...
            'cbh', cbh_l, ...
            'xfov', fov);
    end