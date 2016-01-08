function y = subsref(d,s)

switch s.type
case '()'
  y = feval(d,s.subs{1});
case '.'
  switch s.subs
  case 'x', y = d.x;
  otherwise
    error(['Reference to non-existent field ' s.subs]);
  end;
otherwise
  error('Invalid reference type');
end;
