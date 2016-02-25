function y = subsref(m,s)

switch s.type
case '()'
  y = feval(m,s.subs{1});
case '.'
  switch s.subs
  case 'x', y = m.x;
  case 'kernel', y = m.kernel;
  case 'c', y = m.c;
  case 'a', y = m.a;
  otherwise
    error(['Reference to non-existent field ' s.subs]);
  end;
otherwise
  error('Invalid reference type');
end;
