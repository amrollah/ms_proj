function [ TL ] = prep_LinkeTurbidity( obj )
%PREP_LINK Summary of this function goes here
tmp = load([obj.conf.datafolder 'LinkeTurbidities.mat']);
LinkeTurbidity = double(tmp.LinkeTurbidity);
[~,month,~]=datevec(obj.ti(1));
% Find the appropriate indices for the given Latitude and Longitude
LatitudeIndex = round(LinearlyScale(obj.calib.model3D.Location.latitude, 90, -90, 1, 2160));
LongitudeIndex = round(LinearlyScale(obj.calib.model3D.Location.longitude, -180, 180, 1, 4320));

% Create the "Lookup3D" function to allow fast vector searches into a 3D
% table. This essentially creates a linear index based on the input indices
% and the size of the array, then indexes into the array.
Lookup3D = @(array,a,b,c) array(((c-1).*numel(array(:,:,1))+(b-1).*numel(array(:,1,1))+(a-1)+1));
L1 = Lookup3D(LinkeTurbidity, LatitudeIndex, LongitudeIndex, month);
TL = double(L1)./20;

function OutputMatrix = LinearlyScale(inputmatrix, inputmin, inputmax, outputmin, outputmax)
% OutputMatrix = LinearlyScale(inputmatrix, inputmin, inputmax, outputmin, outputmax)
% Linearly scales the inputmatrix. Maps all values from inputmin to
% outputmin, and from inputmax to outputmax. Linear mapping from one point
% to the other.
    inputrange = inputmax-inputmin;
    outputrange  = outputmax-outputmin;
    
    OutputMatrix = (inputmatrix-inputmin)*outputrange/inputrange + outputmin;

end
end

