function pointReal = GetImage2Real(m,calibModel,ecalibModel,cbh)
% here we project a set of points given in m onto a plane that is cbh meters away from the camera the used convention in opencv is [column,row] kind of
% represantation for pixel points, here however we change that to [row,column]
pointTemp = cam2world([m(:,2),m(:,1)]',calibModel); % we get the corresponding direction on a unit-sphere at first.
pointReal = zeros(size(pointTemp,2),2); % we crate the array that will contain the projected positions
if(size(ecalibModel.R,1)==3)
    pointTemp = ecalibModel.R*pointTemp;
    pointReal = (pointTemp(1:2,:).*repmat((-cbh./pointTemp(3,:)),2,1))';
else
    pointReal = ((ecalibModel.R*pointTemp(1:2,:)).*repmat((-cbh./pointTemp(3,:)),2,1))'; % here we map the point on the unit-sphere onto the plane
end
%for i = 1:size(pointTemp,2)
%    % this part is a bit messy !!!! take care!!!/
%    pointReal(i,:) = (ecalibModel.R*pointTemp(1:2,i))*(-cbh/pointTemp(3,i)); % here we map the point on the unit-sphere onto the plane
%end
%pointReal = mat2cell(pointReal,size(pointReal,1),size(pointReal,2));
end

% % here we project a set of points given in m onto a plane that is cbh meters away from the camera the used convention in opencv is [column,row] kind of
% % represantation for pixel points, here however we change that to [row,column]
% pointTemp = cam2world([m(:,2),m(:,1)]',calibModel); % we get the corresponding direction on a unit-sphere at first.
% pointReal = zeros(size(pointTemp,2),2); % we create the array that will contain the projected positions
% for i = 1:size(pointTemp,2)
%     % this part is a bit messy !!!! take care!!!/
%     pointReal(i,:) = (ecalibModel.R*pointTemp(1:2,i))*(cbh/abs(pointTemp(3,i)))/1; % here we map the point on the unit-sphere onto the plane
%     %pointReal(i,:) = pointTemp(1:2,i)*(cbh/abs(pointTemp(3,i)))/1;
%     %pointReal(i,:) = (ecalibModel.R*pointTemp(1:2,i))*cbh;
% end
% %pointReal = mat2cell(pointReal,size(pointReal,1),size(pointReal,2));