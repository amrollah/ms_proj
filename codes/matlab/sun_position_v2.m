function [ pos, im_pos,Q] = sun_position_v2(obj,j,R)
%SUN_POSITION Summary of this function goes here
    
    [Q, sun] = obj.sunpos_realworld_v2(j);
    %disp(sun.zenith);
    pos = obj.sun_pos([2 3],j);
    %pos = obj.camworld2im(R'*Q);
    %pos = sun_pos_adjuster( pos, pos_m ); 
    im_pos = obj.detect_saturated_sun(j);
    %rmse = sqrt(sum((im_pos-pos).^2));
%     if rmse<obj.conf.sun_detection_thresold
%         pos = im_pos;
%     end

function [ pos ] = sun_pos_adjuster( pos, pos_m )
    %SUN_POS_ADJUSTER Summary of this function goes here
    sign_dist=sign(pos_m-pos);
    dist = sign_dist(1)*pdist([pos,pos_m]');
     % these numbers are set emprically based on RMSE observations
    f1 = abs(dist)-400;
    shift = ((0.00006*(f1.^2)-6)*(1+sign(f1))/2) + ((0.00002*(f1.^2))*(1-sign(f1))/2);
    pos = pos+[1.9*sign_dist(1)*shift; 1.25*sign_dist(1)*shift];
end
end

