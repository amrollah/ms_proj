function y = ftanh(thres,conf,v,have_v)
if any(size(v)~=size(thres.j0))
  if numel(v) == 1, v = zeros(size(thres.j0))+v; else
    x = v; v = have_v; v(logical(have_v)) = x; v(~have_v) = NaN;
  end
end
has_tanh = squeeze(all(~isnan(thres.bb),1));
ally = sum(thres.bb(:,has_tanh).*...
  tanh(bsxfun(@minus,v(has_tanh)',thres.zz')/conf.r2b_resolution));
y = sum(ally);
