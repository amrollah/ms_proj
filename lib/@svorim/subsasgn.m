function m = subsasgn(m,s,x)

switch s.type
case '.'
  switch s.subs
%   case 'lc'
%     if numel(x)~=2, error('lc must be a vector with 2 elements'); end
%     m.lc(1:2) = x;
  otherwise
    error(['Cannot assign to field ' s.subs]);
  end;
otherwise
  error('Invalid reference type');
end;
