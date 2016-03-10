function [ normalized ] = normalize( array, x, y )
%NORMALIZE
     % Normalize to [0, 1]:
     array = double(array);
     mn = min(array,[],2);
     mx = max(array,[],2);
     range = mx - mn;
     array = (array - repmat(mn,1,size(array,2))) ./ repmat(range,1,size(array,2));

     % Then scale to [x,y]:
     range2 = y - x;
     normalized = (array.*range2) + x;
end

