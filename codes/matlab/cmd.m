j = 6300;
j=4800;
obj.loadframe(j)
figure;
obj.showframe(j);
%obj.showThres()
%figure(1);imagesc(obj.xcur.r2b-obj.xcur.thres.vmfi)
%colorbar vert
f=obj.xcur.r2b>=obj.xcur.thres.vmfi;
g=imresize(f,obj.sz);
mask = obj.skymask_wo_sun(j);
mask(isnan(mask)) = false;
g(~mask) = 0;
P = obj.plant_projection_on_image(j);
[x,y]=find(g>-1);
figure(110);imshow(g)
g2 = inpolygon(x,y,P(1,:),P(2,:));
g2 = reshape(g2,size(g));
shadow = g2&g;
figure(120);imshow(g2)
figure(130);imshow(shadow)
fprintf('Shadow area percent: %.0f  \n', 100*sum(shadow)/sum(g2));


for r_scale=1:6
    gr= g;
    obj.conf.rsun_px = obj.conf.rsun_px.*2;
    mask = obj.skymask_wo_sun(j);
    mask(isnan(mask)) = false;
    gr(mask) = 1;
    gr(~mask) = 0;
    figure;
    imshow(gr);
    pause(1);
end