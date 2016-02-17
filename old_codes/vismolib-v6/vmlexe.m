function vmlexe(fname)
fid = fopen(fname);
l = {};
while 1
  strday = fgetl(fid);
  if ~ischar(strday), break; end
  if strday(1)=='2', l = [l {strday}]; end
end
fclose(fid);

load oldpred
vmlconfig_cavriglia;
for i=1:length(l)
  strday = l{i};
  disp([num2str(i) ' / ' num2str(length(l)) ': loading day ' strday]);
  s = vmlSeq(strday,[],VMLCONF);
  s.process(oldpred);
  calc = s.calc;
  save(['calc' strday],'calc');
  oldpred = s.curpred;
end

