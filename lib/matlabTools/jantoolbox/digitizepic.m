function x = digitizepic(fname,xlb,xub,ylb,yub)
if nargin<3
  xlb = 0;
  xub = 1;
end
if nargin<5
  ylb = 0;
  yub = 1;
end
pic = imread(fname);
figure(9746);clf;
image(pic);
nx = size(pic,2);
ny = size(pic,1);
x = ginput;
if ~isempty(x)
  x = [(x(:,1)-1)/(nx-1) (ny-x(:,2))/(ny-1)];
  x(any(x>1.05|x<-.05,2),:) = NaN;
  x = [xlb+(xub-xlb)*x(:,1) ylb+(yub-ylb)*x(:,2)];
end
