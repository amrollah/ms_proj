% mind=min(min(obj.curseg.d2sun));
% maxd=max(max(obj.curseg.d2sun));
% obj.curseg.d2sun2=255*(obj.curseg.d2sun-mind)/(maxd-mind);
predictTime=5; % 1 minute
cloudH=5;

jstart=4340; jend=4370;
jstart2=4324; jend2=4354;

figure; plot(s.calc.avgcc(cloudH,(jstart:jend)+15),s.calc.irr_pred(predictTime,jstart:jend)./s.data.Irr(jstart:jend,2)','.-')
% figure; plot(s.calc.plant_avgcc2p(5,jstart:jend),s.calc.irr_pred(5,jstart:jend)./s.data.Irr(jstart:jend,2)','.-')
title('cloud coverage vs irradiation ratio');
xlim=([0,1]);
ylim([.2,1.2]);
xlabel('Cloud coverage');
ylabel('Irradiation ratio');


figure; h(1)=plot(jstart2:jend2,s.calc.avgcc(cloudH,(jstart2:jend2)),'b.-');
hold on; h(2)=plot(jstart2:jend2,s.calc.irr_pred(predictTime,jstart:jend)./s.data.Irr(jstart:jend,2)','r.-');
legend([h(1),h(2)],{'Cloud coverage'; 'Irradiance ratio'});
