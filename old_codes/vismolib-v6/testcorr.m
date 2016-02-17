function testcorr(obj,j)
obj.loadframe(j);
aa = (1:360)/180*pi;
n = 3;
rr = 21-n:150+n;
prgb=vmlImage2Polar(obj.oi,obj.calc.seg.yx_sun_detected(:,j),aa,rr,double(obj.curseg.x0));

figure(44);
for i=1:3, subplot(3,1,i);imagesc(prgb(:,:,i));colorbar vert;end

pr2b = vmlRed2Blue(prgb);
figure(4);
subplot(2,2,1);image(uint8(prgb(:,1+n:end-n,:)));
subplot(2,2,2);imagesc(max(-2,min(2,(pr2b(:,1+n:end-n)))));colorbar vert

[~,jj] = sort(mean(pr2b(:,1+n:end-n),2));
corr = zeros(2,130); corr0 = [0 0];
for i=1:2
  c = filter(ones(1,1+2*n)/(1+2*n),1,mean(prgb(jj(1:60),:,i)));
  c = c(1+2*n:end);
  corr0(i) = max(c(end-19:end));
  corr(i,:) = max(c,corr0(i));
end

prgb = prgb(:,1+n:end-n,:);
corrf = max(0,min(1,bsxfun(@rdivide,prgb(:,:,2)-corr0(2),max(1,corr(2,:)-corr0(2)))));
prgb(:,:,1) = max(0,prgb(:,:,1)-bsxfun(@times,corrf,corr(1,:)-corr0(1)));
pr2b = vmlRed2Blue(prgb);

figure(4);
subplot(2,2,3);image(uint8(prgb));
subplot(2,2,4);imagesc(max(-2,min(2,(pr2b))));colorbar vert
