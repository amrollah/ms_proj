function [xa,ya] = importpat(file)
f = fopen(file,'rt');
xa = [];
ya = [];
v = [];
npat = [];
nin = [];
nout = [];
n = 0;
state = 0;
while 1
  s = fgets(f);
  if ~ischar(s), break; end;
  if (isempty(npat) | isempty(nin) | isempty(nout))
    if ~isempty(findstr(s,'No. of patterns'))
      npat = str2num(fliplr(strtok(fliplr(s))));
    elseif ~isempty(findstr(s,'No. of input units'))
      nin = str2num(fliplr(strtok(fliplr(s))));
    elseif ~isempty(findstr(s,'No. of output units'))
      nout = str2num(fliplr(strtok(fliplr(s))));
    end; 
  elseif s(1)~='#'
    v = [v str2num(s)];
    if state==0 & length(v)>=nin
      n = n+1;
      xa(:,n) = v(1:nin)';
      state = 1;
      v = [];
    elseif state==1 & length(v)>=nout
      ya(:,n) = v(1:nout)';
      state = 0;
      v = [];
    end;
  end;
end;
if (isempty(npat) | isempty(nin) | isempty(nout))
  warning('Incomplete header!');
elseif state==1
  warning('Less y data than x data!');
elseif n<npat
  warning('Less patterns than specified!');
elseif n>npat
  warning('More patterns than specified!');
end;
