function jmesh(f)
[xx1,xx2] = ndgrid(0:0.01:1,0:0.01:1); xx = [xx1(:) xx2(:)]';
mesh(xx1,xx2,reshape(feval(f,xx),size(xx1)));