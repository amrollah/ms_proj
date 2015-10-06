function pointReal = GetImage2RealTest(m,calibModel,ecalibModel,cbh,theta2,theta3)
% here we project a set of points given in m onto a plane that is cbh meters away from the camera the used convention in opencv is [column,row] kind of
% represantation for pixel points, here however we change that to [row,column]

R1 = [ecalibModel.R,[0;0];0,0,1];
R2 = [1,0,0;0,sin((theta2+90)*pi/180),sin(theta2*pi/180);0,-cos((theta2+90)*pi/180),-cos(theta2*pi/180)];
R3 = [cos(theta3*pi/180),0,cos((theta3+90)*pi/180);0,1,0;sin(theta3*pi/180),0,sin((theta3+90)*pi/180)];

pointTemp = cam2world([m(:,2),m(:,1)]',calibModel); % we get the corresponding direction on a unit-sphere at first.
pointReal = zeros(size(pointTemp,2),2); % we create the array that will contain the projected positions

for i = 1:size(pointTemp,2)
    % this part is a bit messy !!!! take care!!!/
    pointTemp(:,i) = R3*R2*R1*pointTemp(:,i);
    pointReal(i,:) = pointTemp(1:2,i)*(cbh/abs(pointTemp(3,i)))/1; % here we map the point on the unit-sphere onto the plane
    %pointReal(i,:) = pointTemp(1:2,i)*(cbh/abs(pointTemp(3,i)))/1;
    %pointReal(i,:) = (ecalibModel.R*pointTemp(1:2,i))*cbh;
end
%pointReal = mat2cell(pointReal,size(pointReal,1),size(pointReal,2));
end

