%%
clear kfImPl kfPrPl kfImPl_temp kfPrPl_temp;

horizon = 25;
prPred = zeros(2,horizon,length(cents));
prPredV = zeros(2,horizon,length(cents));
prP = zeros(1,length(cents));
imPred = zeros(2,horizon,length(cents));
imPredV = zeros(2,horizon,length(cents));
imP = zeros(1,length(cents));
imObs = zeros(2,horizon,length(cents));
prObs = zeros(2,horizon,length(cents));
for j = 1:size(cents,1)
    if(j==1)
        %SimpleKF(init_x,sigma_p,sigma_q,sigma_r,type,m,time)
        %kfImPl = SimpleKF_v2(cents(j,:),0.6,[0.01,0.01],1e1,'speed',1,1);
        kfImPl = SimpleKF_v2(cents(j,:),0.6,[0.01,0.01,0.01],1e1,'acceleration',1,1);
        %kfPrPl = SimpleKF_v2(centsProj(j,:),0.9,[0.01,0.15],1e1,'speed',cbh,1);
        kfPrPl = SimpleKF_v2(centsProj(j,:),100,10*[0.001,1,5],1e2,'acceleration',cbh,1);
        %kfPrPl = SimpleKF_v2(centsProj(j,:),1e1,[1e1,1e1],2e2,'speed',cbh,1);
    end
    [kfImPl,~,~] = kfImPl.Estimate();
    [kfPrPl,~,~] = kfPrPl.Estimate();
    kfImPl = kfImPl.Update(cents(j,:)',1);
    kfPrPl = kfPrPl.Update(centsProj(j,:)',cbh);
    prP(j) = trace(kfPrPl.P);
    imP(j) = trace(kfImPl.P);
    kfImPl_temp = kfImPl;
    kfPrPl_temp = kfPrPl;
    for k = 1:horizon
        [kfImPl_temp,imPred(:,k,j),~] = kfImPl_temp.Estimate();
        imPredV(:,k,j) = kfImPl_temp.State([2,5]);
        [kfPrPl_temp,prPred(:,k,j),~] = kfPrPl_temp.Estimate();
        %prPredV(:,k,j) = kfPrPl_temp.State(2:2:4);
        prPredV(:,k,j) = kfPrPl_temp.State([2,5]);
        if(j+k <= size(cents,1))
            imObs(:,k,j) = cents(j+k,:)';
            prObs(:,k,j) = centsProj(j+k,:)';
        end
    end
end
part =1;
xobs = prObs(1,part,:); yobs = prObs(2,part,:); xpred = prPred(1,part,:); ypred = prPred(2,part,:);
xdpred = prPredV(1,part,:); ydpred = prPredV(2,part,:);
%figure, plot(xobs(:),yobs(:),xpred(:),ypred(:));
figure, subplot(4,1,1);
plot([1:length(xobs(:))],xobs(:),[1:length(xpred(:))],xpred(:));
subplot(4,1,2);
plot([1:length(xdpred(:))],xdpred(:));
subplot(4,1,3);
plot([1:length(yobs(:))],yobs(:),[1:length(ypred(:))],ypred(:));
subplot(4,1,4);
plot([1:length(ydpred(:))],ydpred(:));
figure, plot(prP);
%subplot(5,1,5);
%a = xobs(2:end)-xobs(1:end-1);
%windowSize = 5;
%a = filter(ones(1,windowSize)/windowSize,1,a(:));
%plot([2:length(xobs(:))],a(:));

%xobs = imObs(1,part,:); yobs = imObs(2,part,:); xpred = imPred(1,part,:); ypred = imPred(2,part,:);
%figure, plot(xobs(:),yobs(:),xpred(:),ypred(:));