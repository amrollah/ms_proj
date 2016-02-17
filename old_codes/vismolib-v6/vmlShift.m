function varargout = vmlShift(iminfo,vyx,varargin)
% [xsh1 [, xsh2, ...]] = vmlShift(iminfo, vyx, x1 [, x2, ..])
% apply lateral movement (shift) to image
% inputs:
%  - iminfo: struct that contains necessary information of the pixel
%            projection to the plane (e.g. <object>.mfi)
%  - vyx: the velocity vector (y component first) in m on the 1m plane
%  - x1, x2,: the image(s) to shift (dimension ny x nx x K, K>=1 is arbitrary)
% output: xsh1, xsh2, ... are the shifted images

jj = round(bsxfun(@minus,bsxfun(@minus,iminfo.ppx,vyx(:)'),...
  iminfo.prg_yx0)/iminfo.dplane)+1;
j1 = (~isnan(jj(:,1)) & jj(:,1)>=1 & jj(:,1)<=iminfo.prg_NN(1) & ...
  jj(:,2)>=1 & jj(:,2)<=iminfo.prg_NN(2));
jsh = ones(length(j1),1);
jsh(j1)=iminfo.prg_j(jj(j1,1)+iminfo.prg_NN(1)*(jj(j1,2)-1));
jsh(jsh<=0) = 1;
varargout = cell(size(varargin));
for i=1:length(varargin)
  sz = [size(varargin{i}) 1];
  varargout{i} = reshape(varargin{i}(bsxfun(@plus,jsh,(0:sz(3)-1)*length(jsh))),sz);
end
