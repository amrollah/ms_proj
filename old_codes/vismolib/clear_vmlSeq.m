disp('clearing all objects of class vmlSeq');
for v = who'
  if strcmp(eval(['class(' v{1} ')']),'vmlSeq');
    for f = eval([v{1} '.open_figs'])
      try
        close(f);
      catch
      end
    end
    disp(['  -> ' v{1}]);
    clear(v{1});
  end
end
disp('clearing class vmlSeq');
clear vmlSeq