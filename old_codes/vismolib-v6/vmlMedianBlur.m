function y = vmlMedianBlur(x,ksize)

%could use median_filter_2d(x,(ksize-1)/2) instead of cv.medianBlur
%more precise (even for 8bit images) but 10x slower

if isa(x,'logical')
  y = logical(cv.medianBlur(uint8(x),'KSize',ksize));  
elseif isa(x,'uint8') || isa(x,'int8') || ...
    (ksize<=5 && (isa(x,'uint16') || isa(x,'int16')))
  y = cv.medianBlur(x,'KSize',ksize);
else
  xmi = min(x(:)); xma = max(x(:));
  if ksize<=5, nbits = 16; x = uint16((x-xmi)/(xma-xmi)*(2^nbits-1));
  else nbits = 8; x = uint8((x-xmi)/(xma-xmi)*(2^nbits-1)); end
  if xma==xmi, y = double(x); else
    y = double(cv.medianBlur(x,'KSize',ksize))/(2^nbits-1)*(xma-xmi)+xmi;
  end
end
