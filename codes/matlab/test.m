close all;

figure;
x1 = -400:20:400;
x = [-600:20:-400,400:20:600];
plot(x,log2(x)/log2(1.75),'-r')
hold on;
plot(x,0.00006*(x.^2)-6,'-g')
hold on;
plot(x1,0.00002*(x1.^2),'-b')
