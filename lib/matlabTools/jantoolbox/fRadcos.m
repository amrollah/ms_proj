function y = radcos(x)
if nargin<1, x = 2; end;
if isempty(x), x = size(x,1); end;
if length(x(:))==1
  if x<0, y = [zeros(x<-2,1); ones(max(-x-3,0)+2,1)]; return; end;
  switch(x)
  case 3
    y=[ 0.1712 0.6079 0.1058;...
        0.7571 0.5422 0.1313;...
        0.3104 0.5217 0.4907];
  case 4
    y=[ 0.4853    0.1385    0.7359    0.0915;...
        0.6480    0.7587    0.3714    0.1385;...
        0.4473    0.2991    0.6362    0.4829;...
        0.4759    0.2771    0.6592    0.3179];
  case 5
    y=[ 0.4116    0.6355    0.1169    0.7853    0.0809;...
        0.6910    0.5110    0.7587    0.3043    0.1433;...
        0.4095    0.5411    0.2931    0.6976    0.4776;...
        0.4315    0.5725    0.2664    0.7168    0.3114;...
        0.4352    0.5777    0.2620    0.7200    0.2836];
  case 6
    y=[ 0.3614    0.5587    0.8131    0.1014    0.0727    0.7102    0.2878;...
        0.7133    0.5907    0.2771    0.7583    0.1468    0.4116    0.0695;...
        0.3860    0.4897    0.7308    0.2894    0.4737    0.6041    0.6086;...
        0.4020    0.5217    0.7501    0.2593    0.3064    0.6302    0.4709;...
        0.4047    0.5271    0.7533    0.2543    0.2786    0.6346    0.4479;...
        0.4058    0.5293    0.7546    0.2522    0.2671    0.6364    0.4385];
  case 7
    y=[ 0.0895    0.0660    0.7525    0.5013    0.6437    0.3243    0.8317    0.2628;...
        0.7577    0.1493    0.3512    0.6367    0.5014    0.7266    0.2644    0.0710;...
        0.2869    0.4708    0.6481    0.4568    0.5465    0.3698    0.7496    0.5939;...
        0.2542    0.3027    0.6703    0.4864    0.5776    0.3807    0.7702    0.4511;...
        0.2487    0.2746    0.6740    0.4913    0.5828    0.3825    0.7736    0.4273;...
        0.2465    0.2631    0.6755    0.4933    0.5850    0.3832    0.7750    0.4175;...
        0.2453    0.2568    0.6764    0.4944    0.5861    0.3837    0.7758    0.4122];
  otherwise
    y = [0.1648 0.6633 0.1117; 0.7562 0.4779 0.1280];
  end;
  return;
end;
[d,n] = size(x);
if d==2
  y = cos(sqrt(sum(x.^2))*9+2)+cos(x(1,:)*11+2)*0.5+sum((x-0.4).^2).^2*15;
else
  rid = ones(1,d-2);
  rin = ones(1,n);
  b = (0:d-3)'*5+1;
  a = 0.6./(1+b);
  y = cos(sqrt(sum(x(1:2,:).^2))*9+2)+...
    cos(x(1,:).*(8+sum(x(3:end,:),1)*8)+2)*0.5+...
    sum((x(1:2,:)-0.4).^2).^2*15+...
    sum((x(3:end,:)-(b(:,rin).*x(rid,:)+1-x(rid+1,:)).*...
    a(:,rin)-0.2).^2,1)*15/sqrt(d-2);
end;
