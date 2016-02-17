cleardays = {'2015_07_19' '2015_08_04' '2015_08_20' '2015_10_24' ...
  '2015_10_31' '2015_11_12' '2016_01_01' '2016_01_17'};
excl = {[] [] [] [550:800] [] [] [400:520 560:610] [400:420 520:540]};
maxpow = nan(1,8); maxpow(3) = 1.07e6; maxpow(end) = 9.6e5;
% xx1={};yy1={};xx2={};yy2={};
% coeff1 = []; coeff2 = [];
% doy = [];
for i=1:8%length(cleardays)
  s = vmlSeq('cavriglia',cleardays{i});
  doy(i) = s.data.day_of_year;
  x = s.getIrrClear; y0 = s.getP;
  if isnan(maxpow(i)), maxp=max(y0); else maxp = maxpow(i); end
  x = x/max(x); y = y0/maxp;
  jnoon = round(mean(find(x==1)));
  tt = 1:length(x);
  t = tt<=jnoon;
  t1 = tt<=jnoon & y>=0.1;
  t1(excl{i}) = false;
  coeff1(:,i) = [1-x(t1)' (1-x(t1)').^2 (1-x(t1)').^3]\(1-y(t1)');
  xx1{i}=x(t1); yy1{i}=y(t1);
  figure(111+i);clf;
  subplot(2,2,1);s.plotpow;
  subplot(2,2,2);
  plot(xx1{i},yy1{i},'b.',...
    x(t),max(0.5*x(t),1-coeff1(1,i)*(1-x(t))-coeff1(2,i)*(1-x(t)).^2-coeff1(3,i)*(1-x(t)).^3),'r');
  title([num2str(i) ' ' cleardays{i}],'interpreter','none');
  subplot(2,2,3);
  plot(tt(t1),y(t1),'b.');
  subplot(2,2,4);
  p2abs(i)=median(y0(t1)./s.getPrelClear(tt(t1)));
  plot(tt(t1),y0(t1)./s.getPrelClear(tt(t1)),'b.');
  drawnow;
end