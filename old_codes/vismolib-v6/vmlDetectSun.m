function [flag,yx_detected,bestalpha_hex,crit] = vmlDetectSun(x0,conf,yx0)
%[crit,flag,yx_detected,bestalpha_hex] = vmlDetectSun(x0,conf,yx0)
% x0 - image
% conf - configuration struct
% yx0 - center for ROI (typically analytical position of the sun)
% flag - classification result
% yx_detected - best guess of true sun position
% bestalpha_hex - best guess for the angle of the sun rays
% crit - vector of all computed criteria

yx = round(yx0);
aa = (1:360)/180*pi; %angles for polar transformations
crit = -ones(1,5); %initialize outputs
yx_detected = [-1 -1];
bestalpha_hex = -1;
if any(yx<1) || any(yx>[size(x0,1);size(x0,2)]) %ROI outside the image
  flag = -1; return;
end
% get ROI for detecting the "black" spot and star
r = conf.r_roi_px+conf.r2_px;
yy = max(1,yx(1)-r):min(size(x0,1),yx(1)+r);
xx = max(1,yx(2)-r):min(size(x0,2),yx(2)+r);
rsq = bsxfun(@plus,(yy-yx(1))'.^2,(xx-yx(2)).^2);
xf = median_filter_2d(x0(yy,xx,:),1,0,1);
my = [0 0 0]; mx = [0 0 0]; diffcrit = [0 0 0];
for i=1:3 %find "black" dot on each channel
  x1 = moving_avg_2d(double(xf(:,:,i)),conf.r0_px,1);
  x2 = x1; x2(rsq>conf.r_roi_px^2) = inf;
  v_in = min(x2(:)); [my1,mx1]=find(x2<=v_in+5);
  my(i) = round(mean(my1)); mx(i) = round(mean(mx1));
  if sqrt(sum(([yy(my(i));xx(mx(i))]-yx0).^2))<=conf.maxerr_px
    rsq=bsxfun(@plus,(yy'-yy(my(i))).^2,(xx-xx(mx(i))).^2);
    v_out = median(x1(rsq>=conf.r1_px^2 & rsq<=conf.r2_px^2));
    diffcrit(i) = v_out-v_in;
  end
end
ii = find(diffcrit>conf.blacklevelthres);
crit(1) = max(diffcrit);
if isempty(ii) %black dot not found
  flag = 1; return;
end

my0 = round(mean(my(ii))); mx0 = round(mean(mx(ii)));
yx_detected = [yy(my0);xx(mx0)]; %position of black dot
iminfo1=[];
iminfo1.sz = size(xf);
[iminfo1.YY,iminfo1.XX] = ndgrid(1:size(xf,1),1:size(xf,2));
rr=conf.r0_px:conf.r2_px;
rsq=bsxfun(@plus,(yy'-yy(my0)).^2,(xx-xx(mx0)).^2);
ishex = [0 0 0]; arel = [0 0 0]; arel0 = [0 0 0];
for i=ii
  x1 = double(xf(:,:,i)); %look for the star patch: it's black ...
  starpatch = rsq<=conf.r1_px^2 & x1<=conf.blacklevelthres;
  starpatch = find_connected(starpatch,my0,mx0); % ... and connected
  [yy1,xx1]=find(starpatch);
  if ~isempty(yy1)
    ch=squeeze(cv.convexHull([xx1(:)-1 yy1(:)-1]));
    if numel(ch)>=6 %at least 3 points in convex hull
      sp_convexhull = ~cv.fillConvexPoly(true(size(x1)),num2cell(ch,2));
      [yy1,xx1]=find(sp_convexhull); %new candidate for center: center of convex hull
      my(i) = round(mean(yy1)); mx(i) = round(mean(xx1));
      if sqrt(sum(([yy(my(i));xx(mx(i))]-yx0).^2))<=conf.maxerr_px
        rsq1 = bsxfun(@plus,(yy'-yy(my(i))).^2,(xx-xx(mx(i))).^2);
        rmax=max(rsq(sp_convexhull)); rmax1=max(rsq1(sp_convexhull));
        if rmax<=rmax1,  %check if new or old center is better
          rsq1=rsq; my(i)=my0; mx(i)=mx0; rmax1=rmax;
        end
        rmax1 = max(conf.r0_px^2,rmax1);
        jclose = rsq1<=conf.r0_px^2;
        sp_convexhull(jclose)=true; starpatch(jclose)=true;
        %find the relative area of the star shape and its convex hull in
        %polar coordinates, this is the main criterion for the star
        [rr1,cc1]=find(vmlImage2Polar(iminfo1,[my(i) mx(i)],aa,rr,double(starpatch)));
        rrsp = rr(accumarray(rr1,cc1,[],@max));
        [rr1,cc1]=find(vmlImage2Polar(iminfo1,[my(i) mx(i)],aa,rr,double(sp_convexhull)));
        rrspf = rr(accumarray(rr1,cc1,[],@max));
        arel0(i) = sum(sp_convexhull(:))/sum(rsq1(:)<=rmax1);
        arel(i) = 1-mean(rrsp.^2)/mean(rrspf.^2);
        %another criterion: number of long args between points defining the convex hull
        [~,jgap]=max(-diff(sort(sqrt(sum((diff([ch;ch(1,:)])).^2,2)),'descend')));
        ishex(i) = (jgap==6);
      end
    end
  end
end
ii = find(arel0>conf.minarel0);
crit(2) = max(arel0);
if ~isempty(ii)
  crit(3) = max(arel(ii));
  crit(4) = max(ishex(ii));
  my0 = round(mean(my(ii))); mx0 = round(mean(mx(ii)));
  yx_detected = [yy(my0);xx(mx0)];
end
star_rmax = min([conf.star_r2_px size(x0,1)-yx_detected(1) ...
  yx_detected(1)-1 size(x0,2)-yx_detected(2) yx_detected(2)-1]);
if star_rmax>conf.star_r2min_px
  yy = yx_detected(1)-star_rmax:yx_detected(1)+star_rmax;
  xx = yx_detected(2)-star_rmax:yx_detected(2)+star_rmax;
  iminfo1 = [];
  iminfo1.sz = star_rmax*2+[1 1];
  [iminfo1.YY,iminfo1.XX] = ndgrid(yy,xx);
  rr = conf.star_r1_px:star_rmax;
  ww = (10.^-((rr-rr(1))'/(rr(end)-rr(1)))); ww = ww/sum(ww);
  hexdrop = [0 0 0];
  for i=3:-1:1
    zz = vmlImage2Polar(iminfo1,yx_detected,aa,rr,double(x0(yy,xx,i)));
    zz = bsxfun(@minus,zz,mean(zz)); zz = bsxfun(@rdivide,zz,std(zz));
    meanzz = mean(zz(:).^2);
    nn = [4 5 8 10 6]; vv=nn;
    for j=1:length(nn)
      n=nn(j);
      zzrep = squeeze(mean(reshape(zz,size(zz,1)/n,n,size(zz,2)),2));
      zz1 = zz-repmat(zzrep,n,1);
      vv(j)=mean(zz1(:).^2)/meanzz;
    end
    b=[nn*0+1;nn]'\vv';
    hexdrop(i) = max(0,b(1)+6*b(2)-vv(end));
  end
  crit(5) = max(hexdrop);
  [~,bestalpha_hex] = max(mean(zzrep(:,floor(size(zzrep,2)/2):end),2));
end
if all(crit<=conf.maxcrit4cov) flag=2;
elseif all(crit>=conf.mincrit4clear) flag=4;
else flag=3;
end
