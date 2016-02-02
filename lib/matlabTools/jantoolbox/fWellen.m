function y = fWellen(x);
ri = ones(1,size(x,2));
c = rem((1:size(x,1)).^3+3,5)'+1;
d = rem((1:size(x,1)).^7+3,11)'/10;
y = sum(sin(c(:,ri).*(x+d(:,ri))),1)+...
  sum(exp(x.*x([2:end 1],:)),1);
