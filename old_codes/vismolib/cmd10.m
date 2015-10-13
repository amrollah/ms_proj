ttr = 0:.01:1;
ttry1 = ttr*0;
ttry0 = ttr*0;
for i=1:size(v1,2)
  jj1 = 1+(ceil(v1(1,i)*100):floor(v1(2,i)*100));
  ttry1(jj1) = ttry1(jj1)+(v1(2,i)-v1(1,i)).^0;
end
for i=1:size(v0,2)
  jj1 = 1+(floor(v0(1,i)*100):ceil(v0(2,i)*100));
  ttry0(jj1) = ttry0(jj1)+1;
end
