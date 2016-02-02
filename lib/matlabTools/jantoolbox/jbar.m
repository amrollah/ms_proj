function jbar(data,b,leg,lpos,bSeperate)

if iscell(data), I=length(data); J=size(data{1},2);
else I=size(data,3); J=size(data,2); end;
m = []; s = []; vmin = []; vmax = []; t=[];
for i=1:I
  if iscell(data), xx = data{i}; end;
  for j=1:J
    if iscell(data), x = xx(:,j); else x = data(:,j,i); end;
    x = x(~isnan(x));
    if isempty(x)
      m(i,j) = NaN;s(i,j) = NaN;t(i,j) = NaN;
      vmin(i,j) = NaN;vmax(i,j) = NaN;
    else
      m(i,j) = mean(x);s(i,j) = std(x);
      t(i,j)=s(i,j)*tinv(0.975,length(x)-1)/sqrt(length(x));
      vmin(i,j) = min(x);vmax(i,j) = max(x);
    end;
  end;
end;

[n,k] = size(m);
if n==1
  if nargin<5 | bSeperate==0
    n = k; k = 1;
    m = m';s = s'; t = t'; vmin = vmin'; vmax = vmax';
  else
    m = [m; zeros(size(m))];
  end;
end;

if nargin<2 | isempty(b), h = bar(m); else h = bar(m,b); end;
if nargin>=3 & ~isempty(leg), 
  if nargin<4, lpos = []; end;
  if length(lpos)<2, legend(leg,lpos); else
    hleg = legend(leg);
    legpos = get(hleg,'position');
    legpos(1:2) = lpos;
    set(hleg,'position',legpos);
  end;
end;
colormap(repmat(1+1/(4*k)-(1:k)'/(2*k),1,3));

for j=1:k
  xd = get(h(j),'xdata');
  for i = 1:n
    x0 = mean(xd(:,i));
    x1 = xd(1,i)-x0;
    xv = x1*0.6;
    xm = x1*0.4;
    hh = line([x0-xv x0+xv x0+xv x0-xv x0-xv], [m(i,j)-s(i,j) ...
        m(i,j)-s(i,j) m(i,j)+s(i,j) m(i,j)+s(i,j) m(i,j)-s(i,j)]);
    set(hh,'Color',[0 0 0]); %,'Linewidth',1
    hh = line([x0-xv x0+xv x0+xv x0-xv x0-xv], [m(i,j)-t(i,j) ...
        m(i,j)-t(i,j) m(i,j)+t(i,j) m(i,j)+t(i,j) m(i,j)-t(i,j)]);
    set(hh,'Color',[0 0 0],'Linewidth',1);
    hh = line([x0-xm x0+xm x0 x0 x0-xm x0+xm], ...
      [vmin(i,j) vmin(i,j) vmin(i,j) vmax(i,j) vmax(i,j) vmax(i,j)]);
    set(hh,'Color',[0.3 0.3 0.3],'Linewidth',1);
  end;
end;