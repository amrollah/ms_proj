function pointImage = GetReal2Image(m,calibModel,ecalibModel,cbh)
% here we project a set of points given in m back onto the image plane. The used convention in opencv is [column,row] kind of
% represantation for pixel points, here however we change that to [row,column]
if(size(ecalibModel.R,1)==3)
    pointTemp = world2cam(ecalibModel.Rinv*[m';-cbh*ones(1,size(m,1))],calibModel)';
else
    pointTemp = world2cam([ecalibModel.R'*m';-cbh*ones(1,size(m,1))],calibModel)'; % we get the corresponding direction on a unit-sphere at first.
end
pointImage = [pointTemp(:,2),pointTemp(:,1)]-1;
end

% % here we project a set of points given in m back onto the image plane. The used convention in opencv is [column,row] kind of
% % represantation for pixel points, here however we change that to [row,column]
% pointTemp = world2cam([ecalibModel.R'*m';-cbh*ones(1,size(m,1))],calibModel)'; % we get the corresponding direction on a unit-sphere at first.
% pointImage = [pointTemp(:,2),pointTemp(:,1)];
