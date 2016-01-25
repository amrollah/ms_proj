% mind=min(min(obj.curseg.d2sun));
% maxd=max(max(obj.curseg.d2sun));
% obj.curseg.d2sun2=255*(obj.curseg.d2sun-mind)/(maxd-mind);

jstart=4324; jend=4374;
figure; plot(s.calc.avgcc(5,jstart:jend),s.calc.irr_pred(5,jstart:jend)./s.data.Irr(jstart:jend,2)','.-')
% figure; plot(s.calc.plant_avgcc2p(5,jstart:jend),s.calc.irr_pred(5,jstart:jend)./s.data.Irr(jstart:jend,2)','.-')
title('cloud coverage vs irradiation ratio');
xlim=([0,1]);
ylim([.2,2]);
xlabel('Cloud coverage');
ylabel('Irradiation ratio');
