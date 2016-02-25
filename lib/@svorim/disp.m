function disp(m)

if isfield(m.kernel,'sig'),
  kernelinfo = [', kernel=Gauss, sigma=' num2str(m.kernel.sig)];
else
  if m.kernel.b>0
    kernelinfo = [', kernel=inhomogeneous polynomial, degree=' num2str(m.kernel.n)];
  else
    kernelinfo = [', kernel=homogeneous polynomial, degree=' num2str(m.kernel.n)];
  end
end
disp(['svorim object in input dimension ' num2str(size(m.x,1)) kernelinfo]);
disp('  readable fields are: x, c, a, kernel');
fprintf('\n');
