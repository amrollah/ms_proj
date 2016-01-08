function show1(s,j)
x = double(s.imread(j));
z = min(1,max(-1,log2(x(:,:,1))-log2(1+x(:,:,3))));
z(~s.sky_mask) = 0;
figure(21);clf;
s.showPredIm(j,120,2,'r'); colorbar vert;
figure(22);clf;
imagesc(z); colorbar vert;
s.showPredIm(j,120,3,'r');
