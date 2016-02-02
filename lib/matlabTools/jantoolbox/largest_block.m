function b = largest_block(s)
if isempty(s)
  b = [];
  return
end
d = (diff(s)==1);
bestlen = 1;
beststart = 1;
len = 1;
start = 1;
for i = 1:length(d)
  if d(i)==1, len = len+1; else
    if len>bestlen
      bestlen = len;
      beststart = start;
    end
    start = i+1;
    len = 1;
  end
end
if len>bestlen
  bestlen = len;
  beststart = start;
end
b = s(beststart:beststart+bestlen-1);
