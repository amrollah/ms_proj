%clear M
outputVideo = VideoWriter('movie1.avi');
outputVideo.FrameRate = 16;
open(outputVideo);
for t=6000:7600
  img = s.imread(t);
  img = vmlColorify(img,~s.oi.sm,2,64);
  img = vmlMedianDownscale(img,3);
  writeVideo(outputVideo,img);
%   figure(1);clf;s.showframe(t);
%   M(t-5999) = getframe;
end
close(outputVideo);
