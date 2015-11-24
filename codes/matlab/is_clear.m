function [ clear ] = is_clear(obj, t)
%IS_CLEAR checks a day log to see if it's a clear sky at time t --date_num (and k
%seconds from past into future)
clear = 0;
[~, t_index] = min(abs(obj.Irr(:,1)-t));
t=obj.Irr(t_index,1);
% sunrise_idx = find(obj.Irr(:,2)>5,1,'first');
% if obj.ClearSkyOrigRef(t_index,2)<100 && abs(obj.Irr(sunrise_idx,1)-t)<(90*60/86400)
%     clear = 1;
%     return;
% else
if obj.ClearSkyOrigRef(t_index,2)>80 && abs(obj.ClearSkyOrigRef(t_index,2) - obj.Irr(t_index,2))/obj.ClearSkyOrigRef(t_index,2) > obj.conf.diff_to_reference_irr_threshold
    % the measured irradiation is so much lower than refernce that we are
    % sure it's cloudy
    return;
end
t_start = t-obj.conf.clear_sky_window/86400;
t_end = t+obj.conf.clear_sky_window/86400;
s_idx = find(obj.Irr(:,1)>t_start,1,'first');
e_idx = t_index;
%e_idx = find(obj.Irr(:,1)<t_end,1,'last');
irr = obj.Irr(s_idx:e_idx,2);
%fprintf('irr length:%d\n',length(irr));
obj.difff(2,t_index) = s_idx;
obj.difff(3,t_index) = e_idx;
if length(irr)<2*obj.conf.irr_comparison_count
    %error('There is no enough irradiation records to determine clear sky.');
    clear = 1;
else
    % I flatten recorded irradiations by deducting(or adding) them from(to)
    % refernce irradiation
    irr = irr - (obj.ClearSkyOrigRef(s_idx:e_idx,2)-obj.ClearSkyOrigRef(t_index));
    irr = sort(irr);
    diff = sum(irr(end-obj.conf.irr_comparison_count+1:end)-wrev(irr(1:obj.conf.irr_comparison_count)));
    if diff < obj.conf.irr_threshold
        clear = 1;
    end
    obj.difff(1,t_index) = diff;
end
end

