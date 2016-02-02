function y = radcos(x)
global BP
y = fRadcos([x;BP(:,ones(1,size(x,2)))]);
