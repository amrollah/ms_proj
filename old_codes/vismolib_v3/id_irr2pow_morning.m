dd = {...
  '2015_07_19' '' ''; ...
  '2015_08_04' '' ''; ...
  '2015_10_30' '8:30' '9:30'; ...
  '2015_11_12' '' ''; ...
  '2015_12_11' '8:30' '10:00'; ...
  };

for i=1:size(dd,1)
  s = vmlSeq(dd{i,1});
  [~,j]=max(s.data.IrrClear(:,2));
  irmax = s.data.IrrClear(j,2);
  t_irmax = s.data.IrrClear(j,1);
  tt = s.data.P(:,1);
  tt(tt>=t_irmax) = [];
  if ~isempty(dd{i,2})
    t1 = datenum(dd{i,1})+(datenum(dd{i,2})-datenum('0:00'));
    t2 = datenum(dd{i,1})+(datenum(dd{i,3})-datenum('0:00'));
    tt(tt>=t1 & tt<=t2) = [];
  end
  %figure(1);clf;plot(tt,s.getIrr(tt),'b',tt,s.getIrrClear(tt)*clear2irr,'r');datetick;grid on;
  figure(1);clf;plot(s.getIrrClear(tt)/irmax,s.getP(tt));grid on;
  1;
end