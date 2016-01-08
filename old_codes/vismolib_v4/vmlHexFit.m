function mse = vmlHexFit(aa,ww,zz)
n = length(aa)/6;
mse = inf;
j = 0;
for i=0:59
  mse1 = ffit(i,i);
  if mse1<mse, mse=mse1; j=i; end
end
l = fminbnd(@(x)ffit(x,j),j-5,j);
r = fminbnd(@(x)ffit(l,x),j,j+5);
mse = ffit(l,r);

  function mse1 = ffit(l,r)
    pp = zeros(n,6); 
    if l==r, pp(l+1,:)=1; else
      l1 = mod(ceil(l),60); r1 = mod(floor(r),60);
      if l1<=r1, pp((l1:r1)+1)=1; else pp(l1+1:end,:)=1; pp(1:r1+1,:)=1; end
      pp(mod(l1-1,60)+1,:) = ceil(l)-l;
      pp(mod(r1+1,60)+1,:) = r-floor(r);
    end
    XX = [aa*0+1 cos(aa*6) sin(aa*6) pp(:)];
    mse1 = sqrt(mean((XX*(XX\zz)-zz).^2))*ww;
  end

end