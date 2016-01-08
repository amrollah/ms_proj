function d = subsasgn(d,s,x)

switch s.type
case '.'
  error('Assignments to eucl2 fields are not allowed');
otherwise
  error('Invalid reference type');
end;
