conf.wsmooth2_m = 0;
conf.wtarget2 = 1;

%tt = 4476:4480;
%tt = [4462 4469 4475 4480 4488];
%tt = 5173:5:5193;

% s = vmlSeq('cavriglia','2015_07_20');
% tt = 5910:1:5930;

%s = vmlSeq('cavriglia','2016_01_16');
%tt = 1880:1:1900;
tt = 2890:1:2893;

tt_sec = (s.data.ti(tt)-s.data.ti(tt(end)))*86400;
sz1 = 60;
seq = []; seq0 = [];
if isempty(seq)
  for i=1:length(tt)
    s.loadframe(tt(i));
    x = s.curseg.cc; x(~s.mfi.sm) = nan;
    x = vmlDownscale(x,sz1);
    x(~isnan(x)) = 2*(x(~isnan(x))>0.5)-1;
    seq0(:,:,i) = x;
    tic;
    x = vmlMakeSurf(x);
    seq(:,:,i) = x;
    disp([i toc]);
  end
end

% tic;
% [seq1,m] = vmlMotSeg2(seq,tt_sec,conf);
% toc

jj = sum(~isnan(seq),3)>=3;
seqa = reshape(seq,numel(seq)/size(seq,3),size(seq,3));
c=vmlPointwiseL1Regression(tt_sec,seqa(jj,:)');
x = nan(size(x));x(jj)=c(1,:);
m = nan(size(x));m(jj)=c(2,:);

figure(6);showmot(m,-tt_sec(1),seq,x);
