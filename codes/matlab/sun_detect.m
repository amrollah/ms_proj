clear all;
close all;
clc;

addpath(genpath('C:\Users\chamsei\Documents\GitHub\ms_proj\old_codes\vismolib'))
vmlconfig_cavriglia;
conf = evalin('base','VMLCONF');
suns = dir('C:\data\cav\sun_detect_test\');
suns = suns(~vertcat(suns.isdir));

nbins=10;
xbins = linspace(0,255,nbins+1);
xbins = xbins+(xbins(2)-xbins(1))/2;
xbins=xbins(1:nbins);
% load('sun_patches.mat');
load('sun_patches_temp.mat');

p_w = 4;
xgrid=1:p_w:40;
ygrid=1:p_w:40;
[X,Y] = meshgrid(xgrid,ygrid);
p_area = (p_w-1)^2;
% 
% for i=1:size(sun_patches,2)/4
%     figure(1);
%     subplot(14,15,i);
%     rgb = sun_patches{1,i};
%     imshow(rgb);
%     sun=imcrop(rgb,sun_patches{2,i});
%     
%     for d=1:3
%         figure(d+4);
%         subplot(14,15,i);    
%         [f,x]=hist(double(reshape(rgb(:,:,d),1,[])),xbins);
%         bar(x,f/sum(f));
%         axis off;
%         xlim([0,260]);
%         ylim([0,1]);
%     end
%     
%     for d=1:3
%         figure(d+1);
%         subplot(14,15,i);    
%         [f,x]=hist(double(reshape(sun(:,:,d),1,[])),xbins);
%         bar(x,f/sum(f));
%         axis off;
%         xlim([0,260]);
%         ylim([0,1]);
%     end
%     
%     tvect = [];
%     for d=1:3
%         figure(d+7);
%         subplot(14,15,i);    
%         [fs,xs]=hist(double(reshape(sun(:,:,d),1,[])),xbins);
%         [fp,xp]=hist(double(reshape(rgb(:,:,d),1,[])),xbins);
%         fval = (fs/sum(fs))./((fp+1e-8)/sum(fp));
%         bar(xs,fval);
%         tvect = [tvect, fval];
%         axis off;
%         xlim([0,260]);
% %         ylim([0,1]);
%     end
%     pause(0.1);
% end


train_count = 0;
test_count = 0;
train_data=[];
train_labels=[];
test_data=[];
test_labels=[];
test_freq=3;
for i=1:size(sun_patches,2)
    disp(i);
    rgb = sun_patches{1,i};
    f=figure(1);
%     set(f,'units','normalized','outerposition',[0 0 1 1]);
    imshow(rgb);
    hold on;
%     rectangle('Position',sun_patches{2,i});
    
    %compute image for automatic detection of the sun
%       x = sum(rgb,3)/(3*255);
%       x = vmlSobelAbs(x);
      xg=rgb2gray(rgb);
      xg=vmlGaussianBlur(xg,conf.gaussblur4sundect_px);
      xg2=xg;
%       [Gmap, Gdir] = imgradient(xg);
%       figure;imshow(xg);
      
%       x=~Gmap&(xg<100);
      gthres=0.4*mean(mean(xg));
      xg(xg2>gthres)=0;
      xg(xg2<=gthres)=255;

        %Boundary Label the Filtered Image
        [L num]=bwlabel(xg);
        STATS=regionprops(L,'all');
        cc=[];
        removed=0;
        %Remove the noisy regions 
        for j=1:num
            dd=STATS(j).Area;
            if (dd < 5 || dd > 400)
                L(L==j)=0;
                removed = removed + 1;
                num=num-1; 
            end
        end
        if removed>0
            [L2 num]=bwlabel(L);
            STATS=regionprops(L2,'all');
        end
        max_area=0;
        for j=1:num
            if STATS(j).Area > max_area
                sun_pos = STATS(j).Centroid;
                max_area=STATS(j).Area;
            end
        end
%         % Trace region boundaries in a binary image.
%         [B,L,N,A] = bwboundaries(L2);
%         %Display results
%         if (length(B)>0)
%         figure;
%         subplot(1,2,1),  imshow(L2);title('BackGround Detected');
%         subplot(1,2,2),  imshow(L2);title('Blob Detected');
%         end
%         hold on;
%         for k=1:length(B),
%             if(~sum(A(k,:)))
%                 boundary = B{k};
%                 plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
%                 for l=find(A(:,k))
%                     boundary = B{l};
%                     plot(boundary(:,2), boundary(:,1), 'g','LineWidth',2);
%                 end
%             end
%         end

%       xg(x~=1)=255;
      imshow(xg);
%       sun_pxs=numel(find(x==1));
%       if sun_pxs<50 || std(find(x==1))>200
      if max_area==0
          disp('No sun is detected.');
          continue;
      end
      hold on;
%       [sortedX,sortingIndices] = sort(xg(:),'ascend');
%       maxValueIndices = sortingIndices(1:floor(0.8*sun_pxs));
%       [p,q] = ind2sub(size(rgb),maxValueIndices);
%       sun_pos = [mean(p);mean(q)];
%       
      plot(sun_pos(1),sun_pos(2),'ro','markersize',10);
      pause(1);
      hold off;
    continue;  
    sun=imcrop(rgb,sun_patches{2,i});
%     sun_area = sun_patches{2,i}(3)*sun_patches{2,i}(4);
    fp=[];
    for d=1:3
        [f,x]=hist(double(reshape(rgb(:,:,d),1,[])),xbins);
        fp(d,:) = f;
    end
    for p=1:numel(X)
%         figure(2);
%         subplot(nbins,nbins,p);

%         last_p=rectangle('Position',[X(p),Y(p),p_w-1,p_w-1]);
        patch = imcrop(rgb,[X(p),Y(p),p_w-1,p_w-1]);
%         imshow(patch);
        tvect = [];
        for d=1:3
            [fs,xs]=hist(double(reshape(patch(:,:,d),1,[])),xbins);
            fval = (fs/sum(fs))./((fp(d,:)+1e-8)/sum(fp(d,:)));
            tvect = [tvect, fval];
        end
        area = rectint(sun_patches{2,i},[X(p),Y(p),p_w-1,p_w-1]);
        
         if mod(i,test_freq)==0
            train_count = train_count+1;
            train_data(train_count,:) = tvect;
            train_patches{train_count}=patch;
         else
             test_count = test_count+1;
            test_data(test_count,:) = tvect;
            test_patches{test_count}=patch;
         end
        if area >= 0.9*p_area % inside of sun area
            if mod(i,test_freq)==0
                train_labels(train_count) = 1;
            else
                test_labels(test_count) = 1;
            end
        else %outside
            if mod(i,test_freq)==0
                train_labels(train_count) = -1;
            else
                test_labels(test_count) = -1;
            end
        end
%         pause(0.1);
%         delete(last_p);
    end
end
train_labels=train_labels';
test_labels=test_labels';
SVMStruct = svmtrain(train_data,train_labels,'kernel_function','rbf');
test_result = svmclassify(SVMStruct,test_data);
acu = 100-100*numel(find(test_labels~=test_result))/numel(test_labels);
save('training_data.mat', 'train_data','train_labels','test_data', 'test_labels', 'SVMStruct');

Sample=test_data;
SampleScaleShift = bsxfun(@plus,Sample,SVMStruct.ScaleData.shift);
Sample = bsxfun(@times,SampleScaleShift,SVMStruct.ScaleData.scaleFactor);
sv = SVMStruct.SupportVectors;
alphaHat = SVMStruct.Alpha;
bias = SVMStruct.Bias;
kfun = SVMStruct.KernelFunction;
kfunargs = SVMStruct.KernelFunctionArgs;
f = kfun(sv,Sample,kfunargs{:})'*alphaHat(:) + bias;
Confidence = f*-1;
