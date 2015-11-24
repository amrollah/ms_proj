function [ irr ] = adjust_reference(obj,irr)
%ADJUST MCCLEAR clear sky model estimated values to compensate for bias
    irr = irr.*obj.conf.irr_scale;
    if obj.conf.adaptive_clearsky_reference
        irr=adaptive_adjust(irr);
    end
    function [ref_irr] = adaptive_adjust(ref_irr)
        for i=1:length(ref_irr)
            if ref_irr(i)<1
                obj.adjuster_ei(i) = 0;
                continue;
            end
            t=obj.Irr(i,1);
            t_past = t-obj.conf.adapt_window_size/86400; % start time of adaptation window
            t_future = t+(obj.conf.adapt_window_size/2)/86400; % end time of adaptation windows
            s_past = find(obj.Irr(:,1)>t_past,1,'first'); % index of start time
            s_future = find(obj.Irr(:,1)<t_future,1,'last'); % index of end time
            ei = sum(((ref_irr(i)-obj.Irr(s_past:i,2)).*mean(obj.is_clear_states(max(1,s_past-2*60):min(i+2*60, length(obj.is_clear_states)))))./log2((i-s_past+2):-1:2)')/(ref_irr(i)*(i-s_past+1)); % calc the average weighted error
%             ei = sign(ei)*min(abs(ei),0.09); % to inforce more smooth result
            ei = min(0,ei);
            obj.adjuster_ei(i) = ei;
            ei = mean(obj.adjuster_ei(1:end));%max(1,i-30*60):i
            ref_irr(i:s_future) = ref_irr(i:s_future).*(1-ei./log2(2:1:(s_future-i+2))'); % adjustment of immediate future
%             ref_irr(i:end) = ref_irr(i:end)*(1-ei); % adjustment of immediate future
%             obj.adjuster_ei(i) = ei;
        end
    end
end

