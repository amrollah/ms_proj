function imagesc_level(x,y,z,level,linestyle,linewidth)
if nargin==1 || isscalar(y),
  if nargin<4, linewidth = 2; else linewidth = level; end
  if nargin<3, linestyle = 'k'; else linestyle = z; end
  if nargin<2, level = 0; else level = y; end
  z = x;
  if iscell(z), zl = z{2}; z = z{1}; else zl = z; end
  x = 1:size(z,2);
  y = 1:size(z,1);
else
  if iscell(z), zl = z{2}; z = z{1}; else zl = z; end
  if nargin<4, level = 0; end
  if nargin<5, linestyle = 'k'; end
  if nargin<6, linewidth = 2; end
end
[yy,xx]=ndgrid(y(:),x(:));
if size(z,3)==3, image(x,y,z); else imagesc(x,y,z); end
hold on;
contour(xx,yy,zl,level+[0 0],linestyle,'linewidth',linewidth);
hold off;
