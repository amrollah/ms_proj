function exportpdf(fname,fheight,yoffset,fwidth,xoffset)
if nargin<2, fheight=[]; end
if nargin<3, yoffset=[]; end
if nargin<4, fwidth=[]; end
if nargin<5, xoffset=[]; end
if isempty(fheight), fheight=1; end
if isempty(yoffset), yoffset = (fheight-1)/2; end
if isempty(fwidth), fwidth=1; end
if isempty(xoffset), xoffset = (fwidth-1)/2; end
ps = get(gcf, 'Position');
ratio = abs((ps(4)-ps(2)) / (ps(3)-ps(1)));
psize = get(gcf, 'papersize');
paperWidth = psize(1);
paperHeight = paperWidth*ratio;

set(gcf, 'papersize', [paperWidth*fwidth paperHeight*fheight]);
set(gcf, 'PaperPosition', [paperWidth*xoffset    paperHeight*yoffset   paperWidth paperHeight]);

print(gcf, '-dpdf', fname);
