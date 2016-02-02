function plotcdf(varargin)
% plotcdf(x1,...)
% plotcdf(x1,...,nbins)
% plot cumulative distribution functions of one or more data sets
% data sets may also be included in cell arrays
% 2015-05-06, ABB.CH-RD.C1, Jan Poland

nbins = 500;
N1 = length(varargin);
if isscalar(varargin{N1})
  nbins = varargin{N1};
  N1 = N1-1;
end
N = N1;
for i=1:N1
  if iscell(varargin{i}), N = N+length(varargin{i})-1; end
end
XX = zeros(nbins,N);
PP = XX;
k = 1;
leg = cell(1,N);
for i=1:N1
  if iscell(varargin{i}), J = length(varargin{i}); else J=1; end
  for j=1:J
    if iscell(varargin{i}), x = varargin{i}{j}; else x = varargin{i}; end
    x = x(~isnan(x));
    if numel(x)<=nbins
      xx = [sort(x(:));nan(nbins-numel(x),1)];
      pp = [((1:numel(x))'-0.5)/numel(x);nan(nbins-numel(x),1)];
    else
      pp = ((1:nbins)'-0.5)/nbins;
      xx = interp1(((1:numel(x))'-0.5)/numel(x),sort(x(:)),pp);
    end
    XX(:,k) = xx;
    PP(:,k) = pp;
    if iscell(varargin{i}), leg{k} = [num2str(i) ' / ' num2str(j)];
    else leg{k} = num2str(i); end
    k = k+1;
  end
end
plot(XX,PP);
legend(leg,'location','best');
