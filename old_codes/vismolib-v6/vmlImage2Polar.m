function ZZ = vmlImage2Polar(iminfo,yxcenter,aa,rr,XX)
ZZ = zeros(length(aa),length(rr),size(XX,3));
for i=1:size(XX,3)
  ZZ(:,:,i)=interp2(iminfo.XX-yxcenter(2),iminfo.YY-yxcenter(1),XX(:,:,i),...
    bsxfun(@times,cos(aa(:)),rr(:)'),bsxfun(@times,sin(aa(:)),rr(:)'));
end
