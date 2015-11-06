function [ pos ] = sun_pos_adjuster( pos, pos_m )
%SUN_POS_ADJUSTER Summary of this function goes here
base = 1.75;
sign_dist=sign(pos_m-pos);
dist = sign_dist(1)*pdist([pos,pos_m]');
f1 = abs(dist-120)-500;  % these numbers are set emprically based on RMSE observations
shift = (log2(f1)/log2(base))*(1+sign(f1))/2;
pos = pos+[1.9*sign_dist(1)*shift; 1.25*sign_dist(1)*shift];

end

