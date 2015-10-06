multiObjectTrackingCV_CamCalibAssistedLabel(obj_r,recorded_tracks,ind);


%% 266,231,224,272,94,185
ind = 231;
[meanCentroids,meanCentroidsProj,centroids,centroidsProj]=multiObjectTrackingCV_CamCalibAssistedLabel(obj_r,recorded_tracks,ind);

smC = sum(meanCentroids,2);
mC = meanCentroids(smC~=0,:);
sC = sum(centroids,2);
C = centroids(sC~=0, :);
smthmC = norm(mean(abs(mC(3:end,:)-2*mC(2:end-1,:)+mC(1:end-2,:))));
smthmCv = norm(std(abs(mC(3:end,:)-2*mC(2:end-1,:)+mC(1:end-2,:))));
smthC = norm(mean(abs(C(3:end,:)-2*C(2:end-1,:)+C(1:end-2,:))));
smthCv = norm(std(abs(C(3:end,:)-2*C(2:end-1,:)+C(1:end-2,:))));

h1 = figure;
hold on;
scatter(mC(:,1),-mC(:,2),50,'MarkerEdgeColor','b','MarkerFaceColor','c');
scatter(C(:,1),-C(:,2),50,'MarkerEdgeColor','r','MarkerFaceColor','y');
title('Object Detection on the Image Plane');
xlabel('Width (pixels)');
ylabel('Height (pixels)');
set(gca,'YTickLabel',num2str(-get(gca,'YTick')'));
legend(strcat('Mean (smth: ',num2str(smthmC,'%.2f'),setstr(177),num2str(smthmCv,'%.2f'),')'),...
    strcat('Bounding Box (smth: ',num2str(smthC,'%.2f'),setstr(177),num2str(smthCv,'%.2f'),')'),...
    'Location','NorthWest');
set(gca,'YTickLabel',num2str(-get(gca,'YTick')'));
hold off;
saveas(h1,strcat('Plots/centroid_sct_',num2str(ind),'.fig'),'fig');
saveas(h1,strcat('Plots/centroid_sct_',num2str(ind),'.jpg'),'jpg');

h2 = figure;
plot(mC(:,1),-mC(:,2),'-bo',C(:,1),-C(:,2),'-rx','LineWidth',1.5,'MarkerSize',6);
title('Object Detection on the Image Plane');
xlabel('Width (pixels)');
ylabel('Height (pixels)');
set(gca,'YTickLabel',num2str(-get(gca,'YTick')'));
legend(strcat('Mean (smth: ',num2str(smthmC,'%.2f'),setstr(177),num2str(smthmCv,'%.2f'),')'),...
    strcat('Bounding Box (smth: ',num2str(smthC,'%.2f'),setstr(177),num2str(smthCv,'%.2f'),')'),...
    'Location','NorthWest');
saveas(h2,strcat('Plots/centroid_plt_',num2str(ind),'.fig'),'fig');
saveas(h2,strcat('Plots/centroid_plt_',num2str(ind),'.jpg'),'jpg');

smCP = sum(meanCentroidsProj,2);
mCP = meanCentroidsProj(smCP~=0,:);
sCP = sum(centroidsProj,2);
CP = centroidsProj(sCP~=0, :);
smthmCP = norm(mean(abs(mCP(3:end,:)-2*mCP(2:end-1,:)+mCP(1:end-2,:))));
smthmCPv = norm(std(abs(mCP(3:end,:)-2*mCP(2:end-1,:)+mCP(1:end-2,:))));
smthCP = norm(mean(abs(CP(3:end,:)-2*CP(2:end-1,:)+CP(1:end-2,:))));
smthCPv = norm(std(abs(CP(3:end,:)-2*CP(2:end-1,:)+CP(1:end-2,:))));

h3 = figure;
hold on;
scatter(mCP(:,1),mCP(:,2),50,'MarkerEdgeColor','b','MarkerFaceColor','c');
scatter(CP(:,1),CP(:,2),50,'MarkerEdgeColor','r','MarkerFaceColor','y');
title('Object Detection in Reality');
xlabel('X (meters)');
ylabel('Y (meters)');
legend(strcat('Mean (smth: ',num2str(smthmCP,'%.2f'),setstr(177),num2str(smthmCPv,'%.2f'),')'),...
    strcat('Bounding Box (smth: ',num2str(smthCP,'%.2f'),setstr(177),num2str(smthCPv,'%.2f'),')'),...
    'Location','North');
hold off;
saveas(h3,strcat('Plots/centroid_proj_sct_',num2str(ind),'.fig'),'fig');
saveas(h3,strcat('Plots/centroid_proj_sct_',num2str(ind),'.jpg'),'jpg');

h4 = figure;
plot(mCP(:,1),mCP(:,2),'-bo',CP(:,1),CP(:,2),'-rx','LineWidth',1.5,'MarkerSize',6);
title('Object Detection in Reality');
xlabel('X (meters)');
ylabel('Y (meters)');
legend(strcat('Mean (smth: ',num2str(smthmCP,'%.2f'),setstr(177),num2str(smthmCPv,'%.2f'),')'),...
    strcat('Bounding Box (smth: ',num2str(smthCP,'%.2f'),setstr(177),num2str(smthCPv,'%.2f'),')'),...
    'Location','North');
saveas(h4,strcat('Plots/centroid_proj_plt_',num2str(ind),'.fig'),'fig');
saveas(h4,strcat('Plots/centroid_proj_plt_',num2str(ind),'.jpg'),'jpg');