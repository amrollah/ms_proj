function show_pow1(s,figno)
figure(figno);clf
tt=s.data.tmin:1/86400:s.data.tmax;
tt1=s.data.tmin:1/86400:(s.data.tmin+s.data.tmax)/2;tt2=tt1(end):1/86400:s.data.tmax;
subplot(1,2,1);
plot(tt,s.getIrr(tt),'b',tt,s.getIrrClear(tt),'r');datetick
grid on;
title(s.uidata.strtitle,'interpreter','none');
subplot(1,2,2);
plot(s.getIrrClear(tt1),s.getP(tt1),'b',s.getIrrClear(tt2),s.getP(tt2),'r');
grid on;