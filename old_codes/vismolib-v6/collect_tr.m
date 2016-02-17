
ksize = s.mfi.ksize*s.mfi.sdf.ksize;
kk = -floor(ksize/2):floor(ksize/2);

p=s.sunpos_im(3882);

[ii,jj] = find(bsxfun(@plus,(s.mfi.sdf.yy'-p(1)).^2,(s.mfi.sdf.xx-p(2)).^2)<=s.conf.rsun_px^2);

i = 1;

round(s.mfi.sdf.xx(jj(i)))+kk