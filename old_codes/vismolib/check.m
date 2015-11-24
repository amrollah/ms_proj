
N = 0;
opts.printlevel = 5;
opts.max_iter = 0;
opts.derivative_test = 'first-order';

Q1 = Q;%+speye(nx)*1e-1;

if 1
  [status,x]=ipoptqp_tanh(x0*0,tanh_spec(1:N,:),...
    Q1,zeros(1,nx),[],[],[],[],....
    zeros(1,nx)-8,zeros(1,nx)+8,opts);
else
  [status,x]=ipoptqp(Q1,zeros(1,nx),[],[],[],[],....
    zeros(1,nx)-8,zeros(1,nx)+8,opts);
end