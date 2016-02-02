function contour0seq(seq)
cols = 'rmgcbk';
clf; hold on;
for i=1:size(seq,3)
  contour(seq(:,:,i),[0 0],cols(i));
end
hold off;
box on;
grid on;
set(gca,'YDir','reverse');
