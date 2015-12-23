function vyx = vmlOpticalFlow(iminfo,conf,x,tt,vyx0,w)
if nargin<5, vyx0 = [0 0]; end
if nargin<6, w = []; end

tt = (tt-tt(end))*86400;
dt = mean(diff(tt));

opts = []; opts.maxFunEval = 40; 
opts.rho_beg=0.05; %opts.xtol_rel = conf.xtol_rel;
ub = conf.vmax*[1 1];
if ~isempty(w), w=repmat(w,[1 1 size(x,3)-1]); end
vyx = bobyqa(@ffit,vyx0,-ub,ub,opts,iminfo,w,x,dt); %conf.vstep*[1 1]
  
  function d = ffit(vyx)
    x0 = x(:,:,2:end);
    xsh = vmlShift(iminfo,vyx*dt,x(:,:,1:end-1));
    j = ~isnan(x0) & ~isnan(xsh);
    if isempty(w)
      d = mean((abs(x0(j)-xsh(j))).^2);
    else
      d = mean((fabs(x0(j)-xsh(j))).^2.*w(j));
    end
  end
end

