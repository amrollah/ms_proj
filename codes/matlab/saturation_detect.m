close all;
clear all;
clc;
training = false;
img_save_path='';
prj_path='';
proj_path;

addpath(prj_path);
s = vml('cavriglia','2015_07_30',[],true);

conf = [];
conf=local_conf(conf);
load('calc\data_clean.mat', 'data');

if training   
    imgs = dir(img_save_path); % get all content of Cavriglia image folder
    imgs = imgs(~vertcat(imgs.isdir)); % filter only folders
    sat_vals = {};
    all_patches = [];

    for i=1:3
        d=data{i};
        if d.sun_flag==4
            continue;
        end
        img_file = [img_save_path d.day '__' num2str(d.j) '.jpeg'];
        I = s.get_image(img_file);
        sun_pos = s.sunpos_im(d.time);
        cx=floor(sun_pos(1)); cy=floor(sun_pos(2)); R=100;
        sun_patch = I(max(cx-R,1):min(cx+R,size(I,1)),max(cy-R,1):min(cy+R,size(I,2)),:); 
        h=figure(1);
    %     maxfig(h,1);
        [x, y, BW, xi, yi] = roipoly(sun_patch);
        rI = reshape(sun_patch,[],3);
        patch = rI(BW,:);
        all_patches = [all_patches; patch];
    %     figure(3); imshow(BW);
        pause();
    end
    all_patches=double(all_patches);
    save('calc\all_patches2.mat','all_patches');
    figure(5);
    subplot(2,2,1);
    rh=histogram(all_patches(:,1),10);
    title('red');
    subplot(2,2,2);
    gh=histogram(all_patches(:,2),10);
    title('green');
    subplot(2,2,3);
    bh=histogram(all_patches(:,3),10);
    title('blue');

else
    load('calc\all_patches2.mat','all_patches');
    [rh,rc]=hist(all_patches(:,1),10);
    [gh,gc]=hist(all_patches(:,2),10);
    [bh,bc]=hist(all_patches(:,3),10);
    rh = rh ./ (size(all_patches,1));
    gh = gh ./ (size(all_patches,1));
    bh = bh ./ (size(all_patches,1));
    p_w = 4;
    R = 60;
    thres=200;
    %evaluation
    for i=1:length(data)
        d=data{i};
        tm = datevec(d.time);
        if d.sun_flag==4 || tm(4)<12 || tm(4)>15
            continue;
        end
        img_file = [img_save_path d.day '__' num2str(d.j) '.jpeg'];
        I = s.get_image(img_file);
        sun_pos = s.sunpos_im(d.time);
        cx=floor(sun_pos(1)); cy=floor(sun_pos(2));
        xgrid = max(cx-R,1):p_w:min(cx+R,size(I,1));
        ygrid = max(cy-R,1):p_w:min(cy+R,size(I,2));
        [X,Y] = meshgrid(xgrid,ygrid);
        area = imcrop(I,[ygrid(1),xgrid(1),xgrid(end)-xgrid(1),xgrid(end)-xgrid(1)]);
        g_area = double(rgb2gray(area));
        mg_area = medfilt2(g_area);
        flat_g_area = g_area(:);
        [areah,areac]=hist(flat_g_area(:),255);
        areah = areah ./ (size(flat_g_area,1));
%         figure(1010); 
%         bar(areah);
      
        sat_border = g_area>200 & g_area<230;
        edge = cv.Canny(g_area,250);
        sat_contour = g_area>220;
        sat_contour = bwareaopen(sat_contour,4);
        dw=zeros(size(g_area));
        u=zeros(size(g_area));
        r=zeros(size(g_area));
        l=zeros(size(g_area));
        inside_sat = false(size(g_area));
        inside_sat2 = false(size(g_area));
        for xi=1:size(g_area,1)
            dw(xi,:) = sum(sat_contour(xi:end,:),1);
            u(xi,:) = sum(sat_contour(1:xi,:),1);
        end
        for yi=1:size(g_area,2)
            r(:,yi) = sum(sat_contour(:,yi:end),2);
            l(:,yi) = sum(sat_contour(:,1:yi),2);
        end
        for xi=1:size(g_area,1)
            for yi=1:size(g_area,2)
                if ~sat_contour(xi,yi) && sum(r(max(1,xi-1):min(xi+1,size(g_area,1)),yi))>0 && sum(l(max(1,xi-1):min(xi+1,size(g_area,1)),yi))>0 && sum(dw(xi,max(yi-1,1):min(yi+1,size(g_area,2))))>0 && sum(u(xi,max(yi-1,1):min(yi+1,size(g_area,2))))>0
                    inside_sat(xi,yi) = true;
                end
                if ~sat_contour(xi,yi) && dw(xi,yi)>0 && u(xi,yi)>0 && r(xi,yi)>0 && l(xi,yi)>0 
                    inside_sat2(xi,yi) = true;
                end
            end
        end
        inside_sat = bwareaopen(inside_sat,10);
        h=figure(123); subplot(2,3,1); imshow(area);
        subplot(2,3,2); imshow(inside_sat);
        subplot(2,3,5); imshow(inside_sat2);
        subplot(2,3,3); imshow(inside_sat & ~inside_sat2); 
        subplot(2,3,6); imshow(sat_contour);
        subplot(2,3,4); imshow(edge);
        maxfig(h,1);
        sat_fact = min(.25,sum(sum(inside_sat))/numel(g_area));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if false
        contours = cv.findContours(sat_contour,'Mode','External');
        curve = cv.approxPolyDP(sat_contour,'Closed', true);
        areas = cell2mat(cellfun(@(x) cv.contourArea(x(:)) ,contours,'UniformOutput',false))';
        % elminate the objects with contour smaller than 350.
        consider = (areas > 150 & areas < 3000); % if the area of an object is less than 300 pixels, we discard that
        % take only the objects that we want to consider
        contours = contours(consider);
        areas = areas(consider);
        im = cv.drawContours(area,contours,'Color',[0 0 0],'Thickness',1);
        figure; imshow(im);
        total_sat_area = sum(areas);
        
%         g_area = imadjust(g_area);
%         I3 = imcomplement(g_area);
%         level = graythresh(mg_area);
%         bw = im2bw(mg_area,level);
        bw = bwareaopen(sat_border,10);
        cc = bwconncomp(bw,8);
        labeled = labelmatrix(cc);
        RGB_label = label2rgb(labeled, @spring, 'c', 'shuffle');
        figure(800); imshow(RGB_label);
        regions_data = regionprops(cc,'basic');
        sun_area = reshape(area,[],3);
        sat_area = false(size(bw));
        for rg=1:length(regions_data)
            if regions_data(rg).Area > 150 && regions_data(rg).Area < 10000
               patch = sun_area(cc.PixelIdxList{rg},:);
               [prh,prc]=hist(patch(:,1),rc);
               [pgh,pgc]=hist(patch(:,2),gc);
               [pbh,pbc]=hist(patch(:,3),bc);
               prh = prh ./ (size(patch,1));
               pgh = pgh ./ (size(patch,1));
               pbh = pbh ./ (size(patch,1));
               distance = pdist2(prh,rh) + pdist2(pbh,bh) + pdist2(pgh,gh) + (pdist2([sun_pos(2),sun_pos(1)],[ygrid(1),xgrid(1)] + regions_data(rg).Centroid)/(.75*R))^2;
               disp(distance);
               if distance < 3.7
                   sat_area(cc.PixelIdxList{rg}) = true;
               end
            end
        end
        figure(500); imshow(sat_area);
        sat_fact = sum(sum(sat_area))/numel(bw);
        
        h=figure(1);
        imshow(area);hold on;
%         rectangle('Position',[ygrid(1),xgrid(1),xgrid(end)-xgrid(1),xgrid(end)-xgrid(1)]);
%         hold on;
%         plot(sun_pos(2),sun_pos(1),'ro','markersize',10);
%         if false
        pch=[];
        pc=1;
        for p=1:numel(X)
            patch = imcrop(I,[Y(p),X(p),p_w-1,p_w-1]);
            r_g = mean(patch(:,1)-patch(:,2));
            g_b = mean(patch(:,2)-patch(:,3));
            [prh,prc]=hist(patch(:,1),rc);
            [pgh,pgc]=hist(patch(:,2),gc);
            [pbh,pbc]=hist(patch(:,3),bc);
            prh = prh ./ (size(patch,1));
            pgh = pgh ./ (size(patch,1));
            pbh = pbh ./ (size(patch,1));
            sr=sort(prh);
            sg=sort(pgh);
            sb=sort(pbh);
%             if max(patch(:,1))>200 || mean(patch(:,3))<100 || mean(patch(:,1))<150 || r_g<1 || r_g>50 || g_b<1 || g_b>50% || sum(sr(end-1:end))<.3 || sum(sg(end-1:end))<.3 || sum(sb(end-1:end))<.3
%                 continue;
%             end
            pch(pc).rh=prh;
            pch(pc).gh=pgh;
            pch(pc).bh=pbh;
            pch(pc).patch=patch;
            pch(pc).pos = [Y(p),X(p)];
            patch=reshape(patch,[],3);
%             distance = compareHists(prh,rh) + compareHists(pbh,bh) + compareHists(pgh,gh) + mean([compareHists(prh,pgh),compareHists(prh,pbh),compareHists(pbh,pgh)]); %(pdist2(prh,rh)+pdist2(pgh,gh)); %+pdist2(pbh,bh)
            distance = pdist2(prh,rh) + pdist2(pbh,bh) + pdist2(pgh,gh);% + mean([pdist2(prh,pgh),pdist2(prh,pbh),pdist2(pbh,pgh)]);
%             distance = pdist2([prh pgh pbh],[rh gh bh]);% + mean([pdist2(prh,pgh),pdist2(prh,pbh),pdist2(pbh,pgh)]);
%             distance = distance + pdist2(sun_pos',[X(p),Y(p)])/(R);
            pch(pc).distance=distance;
            pc=pc+1;
%             if distance<thres
%                 rectangle('Position',[Y(p),X(p),p_w-1,p_w-1]);
%             end
        end
        [~,sx]=sort([pch.distance],'ascend');
        pchs=pch(sx);
        figure(21);
        ncl=20;
        subplot(ncl,4,1);
        bar(rh);
        title('R');
        subplot(ncl,4,2);
        bar(gh);
        title('G');
        subplot(ncl,4,3);
        bar(bh);
        title('B');
        for p=2:min(ncl,length(pchs))
            pt=pchs(p);
            subplot(ncl,4,4*(p-1)+1);
            bar(pt.rh);
            subplot(ncl,4,4*(p-1)+2);
            bar(pt.gh);
            subplot(ncl,4,4*(p-1)+3);
            bar(pt.bh);
            subplot(ncl,4,4*(p-1)+4);
            imshow(pt.patch);
            title(pt.distance);
        end
        h=figure(100);
        imshow(I);hold on;
        rectangle('Position',[ygrid(1),xgrid(1),xgrid(end)-xgrid(1),xgrid(end)-xgrid(1)]);
        hold on;
        for p=1:min(3*ncl,length(pchs))
%             if pchs(p).distance > 2
%                 break;
%             end
            rectangle('Position',[pchs(p).pos(1),pchs(p).pos(2),p_w-1,p_w-1]);
        end
        end
        
        pause(1);
    end
end