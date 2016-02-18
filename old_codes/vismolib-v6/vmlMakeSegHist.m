function thres = vmlMakeSegHist(thres,conf,ix,iy)

sz = size(thres.r2bq);
jjy = thres.yygb(iy)+1:thres.yygb(iy+1);
jjx = thres.xxgb(ix)+1:thres.xxgb(ix+1);
jjy1 = max(1,thres.yygb(iy)+1-conf.rneigh):min(sz(1),thres.yygb(iy+1)+conf.rneigh);
jjx1 = max(1,thres.xxgb(ix)+1-conf.rneigh):min(sz(2),thres.xxgb(ix+1)+conf.rneigh);
xlim = [jjx(1)-thres.xxgb(ix) length(jjx1)-(jjx1(end)-jjx(end))];
ylim = [jjy(1)-thres.yygb(iy) length(jjy1)-(jjy1(end)-jjy(end))];

norm_hist = diff(xlim)*diff(ylim);
dzzh = thres.zzh(2)-thres.zzh(1);

N = length(thres.zz);
X = [ones(size(thres.X,1),1) thres.X];
spcrit = 1;

q0 = thres.r2bq(jjy,jjx); q0 = q0(:);
NN = length(thres.zzh);

%hh0 = hist(q0(q0>0),1:NN)'*conf.scale_hist/norm_hist;
hh0 = arrayfun(@(x)sum(q0==x),1:NN)'*conf.scale_hist/norm_hist;
w_h0 = wlsq(X,hh0,[],[0 conf.alpha*ones(1,N)],[],[],-X,hh0*0);
yy0 = max(0,X*w_h0);
j0 = otsu_local_min(yy0);
if min(sum(yy0(1:j0-1)),sum(yy0(j0+1:end)))>sum(yy0)*0.1
  spcrit = spcrit*single_peak_crit(yy0,j0);
end

pre_cc1 = thres.pre_cc(jjy,jjx);
pre_cc1 = pre_cc1(:);

hh = (arrayfun(@(x)sum(1-abs(pre_cc1(q0==x))), 1:NN)+...
  cumsum(arrayfun(@(x)sum(pre_cc1(q0==x & pre_cc1>0)), 1:NN))+...
  cumsum(arrayfun(@(x)sum(-pre_cc1(q0==x & pre_cc1<0)), 1:NN),2,'reverse'))'*...
  conf.scale_hist/norm_hist;
w_h = wlsq(X,hh,[],[0 conf.alpha*ones(1,N)],[],[],-X,hh*0);

segcrit = vmlCalcSegCrit(thres.r2bq(jjy1,jjx1),NN,...
  xlim,ylim,conf.rneigh,conf.halfgap)'/thres.norm_segcrit;
w_c = wlsq(X,segcrit,[],[0 conf.alpha*ones(1,N)],[],[],-X,hh*0);
maxc = max(0,max(X*w_c));

w = w_h - w_c;
yy = X*w;
[yymin,j0] = min(yy);
j1 = otsu_local_min(yy);
spcrit1 = single_peak_crit(yy,j0);
if min(sum(yy(1:j1-1)),sum(yy(j1+1:end)))>sum(yy)*0.1
  spcrit1 = min(spcrit1,single_peak_crit(yy,j1));
end
spcrit = spcrit*spcrit1;
if yymin<-0.01 && spcrit1==0
  thres.j0(iy,ix)=j0;
  [~,j1] = max(yy(1:j0-1)); [~,j2] = max(yy(j0+1:end)); j2 = j0+j2;
  hh1 = hh*0;
  hh1(1:j1) = yy(j1)-yy(1:j1)+(j1-(1:j1)')*dzzh*0.1;
  hh1(j2:end) = yy(j2)-yy(j2:end)+((j2:length(hh1))'-j2)*dzzh*0.1;
  w = w + wlsq(X,hh1,[],[0 conf.alpha*ones(1,N)],[],[],-X,hh*0);
  thres.v0lo(iy,ix) = thres.zzh(j0);
  thres.v0hi(iy,ix) = thres.zzh(j0);
end

ryy = (max(yy)-yymin)*0.5;
thres.hist_left_open_crit(iy,ix) = 1-min(1,(yy(1)-yymin)/ryy);
thres.hist_right_open_crit(iy,ix) = 1-min(1,(yy(end)-yymin)/ryy);

thres.plot_data{iy,ix} = {hh0 segcrit w_h w_c w};
thres.ww(:,iy,ix) = w(2:end);
thres.maxsegcrit(iy,ix)= maxc;

q0 = double(q0(q0>0)); m = mean(q0); s = std(q0); s1 = max(s,1);
cdfh = (cumsum(hh0)/sum(hh0))*2-1;
jj = max(round(m-5*s1),1):min(round(m+5*s1),NN);

x = bobyqa(@(x)mean((cdfh(jj)-tanh((jj'-x(1))*x(2))).^2),[m 1/s1],[m-2*s1 1/(2*s1)],[m+2*s1 2/s1]);
mse = sqrt(mean((cdfh(jj)-tanh((jj'-x(1))*x(2))).^2));
a = x(2)/dzzh;
thres.p_tanh_cumsum(iy,ix,:) = [a -a*(thres.zzh(1)+(x(1)-1)*dzzh) length(q0)/norm_hist];
mse = (mse-conf.single_peak_mse_range(1))/diff(conf.single_peak_mse_range);
thres.single_peak_crit(iy,ix) = spcrit * (1-max(0,min(1,mse))^2);

end

function spc1 = single_peak_crit(yy,j0)
if nargin<2, j0=[]; end
if isempty(j0), [~,j0] = min(yy); end
if j0==1 || j0==length(yy), spc1=1;
else
  spc1 = 1-max(0,min(1,(min(max(yy(1:j0-1)),max(yy(j0+1:end)))-yy(j0))/0.1));
end
end

