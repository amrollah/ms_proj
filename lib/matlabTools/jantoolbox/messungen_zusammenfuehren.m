function messungen_zusammenfuehren(n,m)

for i=1:n
  load(['m' num2str(i)]);
  T = messung_T;
  t = messung_t;
  M = messung_M;
  FitStat = messung_FitStat;
  xopt = messung_xopt;
  count = messung_i;
  for k=1:m-1
    load(['m' num2str(k*n+i)]);
    T = [T messung_T];
    t = [t messung_t];
    M = cat(3,M,messung_M);
    FitStat = cat(5,FitStat,messung_FitStat);
    xopt = [xopt messung_xopt];
    count = count+messung_i;
  end;
  messung_T = T;
  messung_t = t;
  messung_M = M;
  messung_FitStat = FitStat;
  messung_xopt = xopt;
  messung_i = count;
  save(['neu/m' num2str(i)],'messung_M','messung_T',...
    'messung_i','messung_t','messung_FitStat','messung_xopt');
end;
