function rpcode(j)
if nargin<1, j=0; end;
spc = char(ones(1,j)*' ');
c = dir;
for i=1:length(c)
  n = c(i).name;
  if c(i).isdir
    if ~strcmp(n,'.') & ~strcmp(n,'..')
      cd(n);
      disp([spc 'cd ' pwd]);
      rpcode(j+1);
      cd ..
    end;
  elseif length(n)>=3 & all(lower(n(end-1:end))=='.m')
    if exist(n)==2, disp([spc n]); pcode(n); end;
  end;
end;
