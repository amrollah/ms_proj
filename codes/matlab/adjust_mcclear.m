function [ out_data ] = adjust_mcclear(data)
%ADJUST MCCLEAR clear sky model estimated values to compensate for bias
out_data = data.*1.075;
end

