function y = reconst(hh)
y = [];
for i=1:size(hh,2)
  y = [y zeros(1,hh(2,i))+hh(1,i)];
end
