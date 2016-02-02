function [haxes, hlines] = multplot(xx,yy,linestyles,titles,varargin)
colors = 'bgrcmky';
if nargin<3, linestyles=''; end
if nargin<4, titles={}; end
if nargin==1 || isempty(xx)
  if nargin==1, yy = xx; end
  if iscell(yy)
    N = length(yy);
    sizes = cellfun(@size,yy,'UniformOutput',0);
    sizes = cat(1,sizes{:});
    xdim = find(all([sizes>1;diff(sizes,1)==0],1),1);
    if isempty(xdim), error('found no dimension to for x'); end
    xx=repmat({1:size(yy{1},xdim)},1,N);
  else
    xx = 1:size(yy,1);
    N = size(yy,2);
  end
elseif iscell(yy)
  N = length(yy);
  if ~iscell(xx), xx=repmat({xx},1,N); end
else
  if size(xx,1)==1, xx=xx'; yy=yy'; end
  if size(yy,1)~=size(xx,1)
    if size(yy,2)==size(xx,1), yy=yy';
    else error('xx and yy do not match in size');
    end
  end
  N = size(yy,2);
end
if ~iscell(linestyles), linestyles = repmat({linestyles},1,N); end
haxes = zeros(1,N); hlines = cell(1,N);
clf;
for i=1:N
  haxes(i)=subplot(N,1,i);
  if iscell(yy)
    if iscell(yy{i})
      hold on;
      hlines1 = cell(1,length(yy{i}));
      for j=1:length(yy{i})
        hlines1{j}=plot(xx{i}{j},yy{i}{j},...
          [colors(1+rem(j-1,length(colors))) linestyles{i}],varargin{:});
      end
      hold off;
      hlines{i}=hlines1;
    else
      hlines{i}=plot(xx{i},yy{i},linestyles{i},varargin{:});
    end
  else
    hlines{i}=plot(xx,yy(:,i),linestyles{i},varargin{:});
  end
  if i<=length(titles), title(titles{i}); end
end
linkaxes(haxes,'x');
