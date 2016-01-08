function disp(d)

disp(['eucl2 object in input dimension ' num2str(length(d.x))]);
if length(d.x)<5,
  disp(['  x = [' num2str(d.x(:)') ']']);
else
  disp(['  x : ' num2str(length(d.x)) '-vector']);
end
fprintf('\n');
