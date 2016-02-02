function thres = vmlMakeSegHist(thres,conf,z,ix,iy)
N = length(thres.zz)-1;
nn=hist(z,thres.zzh)';
nncdf = cumsum(nn)/sum(nn);
jlo = find(nncdf>=conf.outside_quantil,1);
jhi = find(nncdf<1-conf.outside_quantil,1,'last')+1;

X = [ones(size(thres.X,1),1) thres.X];
w0 = wlsq(X,nn,[],[0 conf.alpha*ones(1,N+1)]);
yy0 = max(0,X*w0);
j0 = otsu_local_min(yy0);
scale = sum(nn(jlo:jhi))/thres.scale/sum(yy0(jlo:jhi)); 
lpen = 0; rpen = 0;
if min(sum(yy0(1:j0-1)),sum(yy0(j0+1:end)))<sum(yy0)*conf.min_peak_mass
  %no two peaks detected
  j0 = NaN;
  %scale = scale*conf.scale_singlepeak;
  if thres.localvar(iy,ix)<conf.locvarpen.thres(1) && ...
      thres.r2brange(iy,ix)<=conf.locvarpen.maxrange4sky
    lpen = conf.locvarpen.wsky;
  elseif thres.localvar(iy,ix)>conf.locvarpen.thres(2)
    rpen = conf.locvarpen.wcloud;
  end
elseif jhi-jlo>=round(conf.r2b_maxrange_singlepeak/conf.r2b_hist_resolution)
  %two peaks and even wide enough: penalize outside range
  lpen=1; rpen=1;
end
w0 = w0*scale; nn = nn*scale; yy0 = yy0*scale;
nn1 = nn;
v_out = max(yy0); %mean(nn(jlo:jhi));
if lpen>0
  jlo = jlo+min(find(diff(yy0(jlo:end))<0,1),find(yy0(jlo:end)>=v_out*lpen*0.99,1))-1;
  jj = 1:min(jlo,find(thres.zzh>=conf.r2b_lpenalty_lim,1))-1;
  nn1(jj) = min(nn1(jj)+v_out*lpen,v_out);
end
if rpen>0
  jhi = jhi-min(find(diff(flipud(yy0(1:jhi)))<0,1),find(flipud(yy0(1:jhi))>=v_out*lpen*0.99,1))+1;
  jj = max(jhi,find(thres.zzh<=conf.r2b_rpenalty_lim,1,'last'))+1:length(nn1);
  nn1(jj) = min(nn1(jj)+v_out*rpen,v_out);
end
jj = thres.zzh<conf.r2b_prior(1); 
nn1(jj) = nn1(jj)+(conf.r2b_prior(1)-thres.zzh(jj))'*conf.penalty_prior;
jj = thres.zzh>conf.r2b_prior(end); 
nn1(jj) = nn1(jj)+(thres.zzh(jj)-conf.r2b_prior(end))'*conf.penalty_prior;
w = wlsq(X,nn1,[],[0 conf.alpha*ones(1,N+1)],[],[],-X,nn*0);
thres.j0(iy,ix) = j0;
thres.histn(:,iy,ix) = nn;
thres.bb0(iy,ix) = w0(1);
thres.ww0(:,iy,ix) = w0(2:end);
thres.bb(iy,ix) = w(1);
thres.ww(:,iy,ix) = w(2:end);
