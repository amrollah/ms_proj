function test1(obj,j)

sat = 0.5;
rfade_px = 100;
rmin_px = 20;

obj.loadframe(j);
yxsun = obj.sunpos_im(j);

xx = 1:obj.oi.sz(2); yy = 1:obj.oi.sz(1);
jy = abs(yy-yxsun(1))<=170;
jx = abs(xx-yxsun(2))<=170;
figure(100);clf;
subplot(1,2,1); image(obj.curseg.x0(jy,jx,:));

jy = abs(obj.mfi.yy-yxsun(1))<=170;
jx = abs(obj.mfi.xx-yxsun(2))<=170;
x = obj.curseg.x(jy,jx,:);
d2sun = sqrt((obj.mfi.YY(jy,jx)-yxsun(1)).^2+(obj.mfi.XX(jy,jx)-yxsun(2)).^2);
have_v = d2sun<obj.conf.rsun_px & d2sun>=rmin_px;

rho = max(-sat,min(sat,obj.curseg.r2b(jy,jx)-obj.curseg.thres.vmfi(jy,jx)));
rho = max(-1,min(1,rho/median(abs(rho(have_v)))));
rho(isnan(obj.curseg.r2b(jy,jx))) = nan;

rho(d2sun<rfade_px) = rho(d2sun<rfade_px).*(d2sun(d2sun<rfade_px)/rfade_px);
rho(d2sun<rmin_px) = -1;
subplot(1,2,2);imagesc(rho);colorbar vert;
cc = vmlChanVeseClose2Sun(have_v,x,rho);

% figure(100);clf;
% subplot(1,3,1);imagesc(x);hold on; contour(r2b-r2bthres,[0 0],'k'); hold off;
% subplot(1,3,2);imagesc(x);hold on; contour(ysvm,[0.5 0.5],'k'); hold off;
% subplot(1,3,3);imagesc0(max(-1,min(1,r2b-r2bthres)));
