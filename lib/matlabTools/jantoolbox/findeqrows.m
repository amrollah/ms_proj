function [xi, yi] = findeqrows(x,y)
xi = []; yi = [];
for i=1:size(x,1)
  k = find(all(x(zeros(1,size(y,1))+i,:)==y,2));
  xi = [xi;i(ones(length(k),1))];
  yi = [yi;k];
end;
