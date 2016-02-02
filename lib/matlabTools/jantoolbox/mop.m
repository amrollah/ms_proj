function C = mop(A,op,B)
d = max(ndims(A),ndims(B));
ii = ones(1,d);
for i=1:d
  if size(A,i)~=size(B,i)
    if size(A,i)==1
      ii1 = ii; ii1(i) = size(B,i); A = repmat(A,ii1);
    elseif size(B,i)==1
      ii1 = ii; ii1(i) = size(A,i); B = repmat(B,ii1);
    else
      error(['cannot align dimension ' num2str(i)]);
    end
  end
end
C = eval(['A' op 'B']);
