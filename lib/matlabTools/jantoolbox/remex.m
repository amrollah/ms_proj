function remex(j)
if nargin<1, j=0; end;
spc = char(ones(1,j)*' ');
c = dir;
for i=1:length(c)
  n = c(i).name;
  if c(i).isdir
    if ~strcmp(n,'.') & ~strcmp(n,'..')
      cd(n);
      disp([spc 'cd ' pwd]);
      remex(j+1);
      cd ..
    end;
  elseif length(n)>=3 & all(lower(n(end-1:end))=='.c')
    if exist(n(1:end-2))==3, disp([spc n]); mex(n); end;
  end;
end;
