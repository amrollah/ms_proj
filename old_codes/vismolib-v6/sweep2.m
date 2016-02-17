ee = [0 10.^[-5:.5:-1.5]];
ttpred = [1 2 3 4 5 7 10 15]*60;
% yp = zeros(length(ee),length(ttpred));
% yirr = yp;
% ypp = zeros(1,length(ttpred)); yirrp = ypp;
s.uidata.printlevel = 0;
tt = 5000:5600;
for k=4 %1:length(ttpred)
  for i=1:length(ee)
    disp([i length(ee) k length(ttpred)]);
    s.conf.pred.eta_p=ee(i);
    s.conf.pred.eta_irr=ee(i);
    s.conf.pred.tpred = ttpred(k);
    s.process([],5000,5600);
    
    ccalc{i,k} = s.calc;
    
    d = (abs(s.calc.pred.p.v(tt)-s.getP(tt,ttpred(k))));
    yp(i,k) = mean(d(~isnan(d)));
    ypp(k) = mean(abs(s.getPersistErrP(tt(~isnan(d)),ttpred(k))));
    
    d = (abs(s.calc.pred.irr.v(tt)-s.getIrr(tt,ttpred(k))));
    yirr(i,k) = mean(d(~isnan(d)));
    yirrp(k) = mean(abs(s.getPersistErrIrr(tt(~isnan(d)),ttpred(k))));
    
    disp([yp(i,k)/ypp(k) yirr(i,k)/yirrp(k)]);
    save sweep2 yp yirr ypp yirrp
  end
end
