function FitnessCurves(yrange, varargin)

if rem(length(varargin),5)~=0
  error('Incomplete function call');
end;

figure(26); 
clf;
hold on;

for k = 1:length(varargin)/5
  n = varargin{k*5-4};
  w = varargin{k*5-3};
  ls = varargin{k*5-2};
  xpos = varargin{k*5-1};
  ms = varargin{k*5};
  
  load(['m' num2str(n)]);
  if w==0
    f = mean(min(messung_FitStat(1,:,:,:,:),[],4),5);
  elseif w>0
    f = max(min(messung_FitStat(1,:,:,:,:),[],4),[],5);
  else
    f = min(min(messung_FitStat(1,:,:,:,:),[],4),[],5);
  end;
  
  t = length(f);
  if isempty(xpos), xpos=rand*0.8+0.1; end;
  
  plot([1:t],f,ls);
  plot(round(xpos*t),f(round(xpos*t)),ms);
end;

hold off;
if ~isempty(yrange) axis([0 t yrange]); end;
