function y = quantil(x,p,dim)

if isempty(x), error('x must not be empty'); end
  
sz = size(x);
if nargin < 3 
    dim = find(sz ~= 1,1);
    if isempty(dim)
        dim = 1; 
    end
    dimArgGiven = false;
else
    nDimsX = ndims(x);
    perm = [dim:max(nDimsX,dim) 1:dim-1];
    x = permute(x,perm);
    if dim > nDimsX
        sz = [sz ones(1,dim-nDimsX)];
    end
    sz = sz(perm);
    dim = 1;
    dimArgGiven = true;
end


nrows = sz(dim);
ncols = prod(sz) ./ nrows;
x = reshape(x, nrows, ncols);

x = sort(x,1);
nonnans = ~isnan(x);

% If there are no NaNs, do all cols at once.
if all(nonnans(:))
  n = sz(dim);
  if isequal(p,50) % make the median fast
    if rem(n,2) % n is odd
      y = x((n+1)/2,:);
    else        % n is even
      y = (x(n/2,:) + x(n/2+1,:))/2;
    end
  else
    q = [0 (0.5:(n-0.5))./n 1]';
    xx = [x(1,:); x(1:n,:); x(n,:)];
    y = zeros(numel(p), ncols, class(x));
    y(:,:) = interp1q(q,xx,p(:));
  end
  
  % If there are NaNs, work on each column separately.
else
  % Get percentiles of the non-NaN values in each column.
  y = nan(numel(p), ncols, class(x));
  for j = 1:ncols
    nj = find(nonnans(:,j),1,'last');
    if nj > 0
      if isequal(p,50) % make the median fast
        if rem(nj,2) % nj is odd
          y(:,j) = x((nj+1)/2,j);
        else         % nj is even
          y(:,j) = (x(nj/2,j) + x(nj/2+1,j))/2;
        end
      else
        q = [0 (0.5:(nj-0.5))./nj 1]';
        xx = [x(1,j); x(1:nj,j); x(nj,j)];
        y(:,j) = interp1q(q,xx,p(:));
      end
    end
  end
end

szout = sz; szout(dim) = numel(p);
y = reshape(y,szout);

if dimArgGiven
     y = ipermute(y,perm);  
end

if ~dimArgGiven && isvector(x)
    y = reshape(y,size(p)); 
end
