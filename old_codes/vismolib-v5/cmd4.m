conf.wsmooth2 = 0.01;
conf.wsmooth2_m = 1;
conf.wtarget2_m = 1;
conf.wlo = .01;
conf.whi = 1;
conf.vsat = 0.5;
conf.ub = 10;


%tt = 4476:4480;
%tt = [4462 4469 4475 4480 4488];
%tt = 5173:5:5193;

%tt = 5910:4:5930;
tt = 5910:10:5930;
%tt = [5910 5918];
tt_sec = (s.data.ti(tt)-s.data.ti(tt(end)))*86400;
scale = 0.2;
seq = [];
if isempty(seq)
  for i=1:length(tt)
    s.loadframe(tt(i));
    
    x = replace_isolated_nan(2*s.curseg.cc-1);
    x(~s.mfi.sm) = nan;
    x = cv.resize(x,scale);
    x(~isnan(x)) = 2*(x(~isnan(x))>0)-1;
    seq(:,:,i) = x;
  end
end

tic;
[seq1,m] = vmlMotSeg3(seq,tt_sec,conf);
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
