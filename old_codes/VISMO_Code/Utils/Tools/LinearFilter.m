%% Generic filter implementation.
% c : length of the filter
% d : dimensions of the filter
% functions for initializing a filter, we can build a separate
% class out of this part...
% B. Zeydan, 27. Feb. 2013
%
%
%
% Turned into a class. 
% Burak Zeydan, 22. Aug. 2014

classdef LinearFilter
    
    properties (GetAccess=public)
        capacity
        data
        cursor
        numel
        value
    end
    
    methods
        function obj = LinearFilter(c,d)
            obj.capacity = c;
            obj.data = zeros(c,d);
            obj.cursor = 1;
            obj.numel = 0;
            obj.value = zeros(1,d);
        end
        
        function obj = Add(obj, newData)
            obj.data(obj.cursor,:) = newData;
            obj.cursor = mod(obj.cursor + 1, obj.capacity)+1;
            if(obj.numel < obj.capacity)
                obj.numel = obj.numel + 1;
            end
            obj.value = sum(obj.data,1)/obj.numel;
        end
    end
    
end

