function M = messungen_sammeln2(file,varname)

if nargin==0 | isempty(file), file = 'm';
elseif isnumeric(file), file = ['m' num2str(file)];end;

M = [];
i = 1;
while exist([file num2str(i) '.mat'])==2
  load([file num2str(i)]);
  x = eval(varname);
  M = [[M;zeros(length(x)-size(M,1),size(M,2))+NaN] ...
      [x(:);zeros(size(M,1)-length(x),1)+NaN]];
  i = i+1;
end;
