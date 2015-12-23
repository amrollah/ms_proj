function d = ffit2(vyx,iminfo,w,x,dt)
    x0 = x(:,:,2:end);
    xsh = vmlShift(iminfo,vyx*dt,x(:,:,1:end-1));
    j = ~isnan(x0) & ~isnan(xsh);
    if isempty(w)
      d = mean((abs(x0(j)-xsh(j))).^2);
    else
      d = mean((fabs(x0(j)-xsh(j))).^2.*w(j));
    end
end

