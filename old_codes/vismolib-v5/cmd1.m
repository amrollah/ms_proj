conf.wsmooth = 20;
conf.wsmooth1_m = 0;
conf.wsmooth2_m = 1;
conf.wzero2_m = 1000;
conf.wtarget1 = 0;
conf.wtarget2 = 1;
conf.band = 1;

%tt = 4476:4480;
%tt = [4462 4469 4475 4480 4488];
%tt = 5173:5:5193;
%tt = 5910:4:5930;
tt = [5910 5918];
tt_sec = (s.data.ti(tt)-s.data.ti(tt(end)))*86400;
scale = 0.3;
for i=1:length(tt)
  s.loadframe(tt(i));
  %x = s.curseg.r2b-s.curseg.thres.vmfi;
  %x(~s.curseg.sm) = NaN;
  x = s.curseg.cc*2-1;
  if i==1
    seq = cv.resize(x,scale);
    seq = repmat(seq,[1 1 length(tt)]);
  else
    seq(:,:,i) = cv.resize(x,scale);
  end
end
tic;
[seq1,m] = vmlMotSeg(seq,tt_sec,conf);
toc

% figure(1); clf;
% for i=1:length(tt)
%   subplot(2,length(tt),i);
%   imagesc(seq(:,:,end));
%   colorbar vert;
%   subplot(2,length(tt),i+length(tt));
%   imagesc(seq1(:,:,end));
%   colorbar vert;
% end

figure(3);imagesc(m);colorbar vert
hold on;
contour(seq1(:,:,1),[0 0],'k');
contour(seq1(:,:,end),[0 0],'r');
contour(seq(:,:,end),[0 0],'r:');
%contour(seq1(:,:,end),[0 0],'g');
contour(seq1(:,:,end)-1*tt_sec(1)*m,[0 0],'g');
hold off;
set(gca,'clim',[min(m(:)) max(m(:))]);
