function H = entropy8bit(x)
f = countentries(x+1,256);
f = (f(f>0))/length(x);
H = sum(f.*log(f));
