function dispstruct(s,pre)
if nargin<2, pre = 's'; end
if isempty(s), return; end
if ~isstruct(s), error('argument s must be a struct'); end
fn = fieldnames(s);
for i=1:length(fn)  
  if isstruct(s.(fn{i})), dispstruct(s.(fn{i}),[pre '.' fn{i}]);
  elseif isnumeric(s.(fn{i})) disp([pre '.' fn{i} ' = ' num2str(s.(fn{i}))]);
  else disp([pre '.' fn{i} ' = ' s.(fn{i})]);
  end
end

    