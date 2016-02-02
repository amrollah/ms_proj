function jj = contiguous_idz(x)
jj = x;
i = 0;
k = 1;
while i<length(x)
  j = i+find(isnan(x(i+1:end)),1);
  if isempty(j), j=length(x)+1; end
  jj(i+1:j-1)=k;
  k=k+1;
  i=j;
end
