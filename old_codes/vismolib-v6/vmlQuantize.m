function Y = vmlQuantize(X,xx)
Y = uint8(max(1,min(length(xx),round((X-xx(1))/(xx(2)-xx(1)))+1)));
Y(isnan(X))=0;
