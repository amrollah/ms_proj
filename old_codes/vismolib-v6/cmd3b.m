
Nx = 50; Ny = 50;
r2 = bsxfun(@plus,((1:Ny)'-(1+Ny)/2).^2,((1:Nx)-(1+Nx)/2).^2);
x0 = 2*(r2<=100)-1;
tic;
x = vmlMakeSurf(x0);
toc
figure(99);mesh(x);
