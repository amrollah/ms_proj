load allP
folders0 = folders;
folders0 = folders0(1:length(allP));
folders = dir(VMLCONF.basefolder);
folders = {folders.name};
folders(end)=[];
folders(cellfun(@length,folders)~=10)=[];
%allP = cell(1,length(folders));
%Pncol = zeros(1,length(folders));
for i=length(folders0)+1:length(folders)
  try
    s = vmlSeq(folders{i});
    Pncol(i) = s.data.P_ncols;
    allP{i} = s.data.P;
  catch
    Pncol(i) = NaN;
    allP{i} = [];
  end
end