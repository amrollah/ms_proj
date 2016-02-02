function deobj(cname)

exclude = {'disp.m', 'display.m', 'fcnchk.m', 'feval.m', ...
    'subsasgn.m', 'subsref.m', 'struct.m'};
deobjdir = 'deobjectized';
privdir = 'private';

cpath = which(cname);
cpath = strrep(cpath,'\','/');
j = findstr(cpath,'/');
cpath = cpath(1:j(end));
dest = cpath(1:j(end-1));
if exist([dest deobjdir])~=7, mkdir(dest,deobjdir); end;
dest = [dest deobjdir '/'];
if exist([dest cname])~=7, mkdir(dest,cname); end;
dest = [dest cname '/'];

if exist([cpath privdir])==7
  disp('Copying private...');
  if exist([dest privdir])~=7, mkdir(dest,privdir); end;
  d = dir([cpath privdir]);
  for i=1:length(d) 
    if d(i).isdir==0
      disp(['  ' d(i).name]);
      copyfile([cpath privdir '/' d(i).name], ...
        [dest privdir  '/' d(i).name]);
    end; 
  end;  
end;

disp('Processing methods...');
meth = {};
d = dir(cpath);
for i=1:length(d)
  if d(i).isdir==0
    bex = 0;
    for k=1:length(exclude)
      if strcmp(d(i).name,exclude{k}), bex=1; end;
    end;
    if ~bex, meth{length(meth)+1} = d(i).name; end;
  end;
end;
for i=1:length(meth)
  m = meth{i};  
  if lower(m(end-1:end))~='.m'
    disp(['  ' m ' (copy)']);
    copyfile([cpath m], [dest cname '_' m]);
  else
    disp(['  ' m ' (parse)']);
    in = fopen([cpath m],'rt');
    out = fopen([dest cname '_' m],'wt');
    while 1
      s = fgets(in);
      if ~ischar(s), break; end;
      if ~isempty(findstr(s,'class(')) & ~isempty(findstr(s,['''' cname '''']))
        s = ['% (class construction!) ' s];
      elseif ~isempty(findstr(s,'isa(')) & ~isempty(findstr(s,['''' cname '''']))
        j = findstr(s,['''' cname '''']);
        for l=1:length(j)
          jj = j(l);
          s(jj+1:jj+length(cname)) = [];
          s = [s(1:jj) 'struct' s(jj+1:end)];
          j(l+1:end) = j(l+1:end)-length(cname)+6;
        end;
      else
        for k=1:length(meth)
          m2 = meth{k};
          if length(m2)>2 & (lower(m2(end-1:end))=='.m' | lower(m2(end-1:end))=='.c')
            j = findstr(s,[m2(1:end-2) '(']);
            for l=1:length(j)
              jj = j(l);
              if jj==1 | (~isletter(s(jj-1)) & s(jj-1)~='_' & ~(s(jj-1)>='0' & s(jj-1)<='9'))
                s = [s(1:jj-1) cname '_' s(jj:end)];
                j(l+1:end) = j(l+1:end)+(length(cname)+1);
              end;
            end;
          end;
        end;
      end;
      fprintf(out,'%s',s);
    end;
    fclose(in);
    fclose(out);
  end;
end;
