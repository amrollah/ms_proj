function y = vmlExtend1Px(x)
y = false(size(x));
y(1:end-1,:) = y(1:end-1,:) | x(2:end,:);
y(2:end,:) = y(2:end,:) | x(1:end-1,:);
y(:,1:end-1) = y(:,1:end-1) | x(:,2:end);
y(:,2:end) = y(:,2:end) | x(:,1:end-1);

