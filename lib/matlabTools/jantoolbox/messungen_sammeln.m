function [M, Mstat, T, opt, x_opt] = messungen_sammeln(file)

if nargin==0 | isempty(file), file = 'm';
elseif isnumeric(file), file = ['m' num2str(file)];end;

M = [];
Mstat = [];
T = [];
opt = Inf;
x_opt = [];
i = 1;
while exist([file num2str(i) '.mat'])==2
  load([file num2str(i)]);
  M = [[M;zeros(length(messung_M)-size(M,1),size(M,2))+NaN] ...
      [messung_M(:);zeros(size(M,1)-length(messung_M),1)+NaN]];
  Mstat=[Mstat [min(messung_M);mean(messung_M);max(messung_M)]];
  T = [[T;zeros(length(messung_T)-size(T,1),size(T,2))+NaN] ...
      [messung_T(:);zeros(size(T,1)-length(messung_T),1)+NaN]];
  if messung_opt<opt
    opt = messung_opt;
    x_opt = messung_xopt;
  end;
  i = i+1;
end;
