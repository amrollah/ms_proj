vmlconfig; s.conf=VMLCONF;
[w,b] = s.trainsvm(120);
[yy,bb]=s.predsvm(120,w,b(1));
yypersist = s.persist_cloudiness(120);
figure(11);subplot(2,1,1);plot(b);subplot(2,1,2);plot(bb)
figure(12);s.plot_class_pred(120,yy);
figure(13);s.plot_class_pred(120,yypersist);