targets = {1 2 3 4 5};
x = {0 1 2 3 4}
net = narxnet(1:2,1:2,10);
view(net)
[Xs,Xi,Ai,Ts,EWs,shift] = preparets(net,x,{},targets)