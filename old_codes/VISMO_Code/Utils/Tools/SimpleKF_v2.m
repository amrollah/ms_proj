classdef SimpleKF_v2
    %SIMPLEKF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess=public)
        P
        sigma_p
        x
        R
        sigma_r
        Q
        sigma_q
        m
        A
        H
        A_discr
        time
        State
        type
        dims
        state_to_x
    end
    
    methods
        function obj = SimpleKF_v2(init_x,sigma_p,sigma_q,sigma_r,type,m,time)
            switch type
                case 'force'
                    obj.dims = 6;
                    obj.type = type;
                    obj.P = sigma_p*eye(obj.dims);
                    obj.sigma_p = sigma_p;
                    obj.State = [init_x(1);0;0;init_x(2);0;0];
                    obj.state_to_x = [1:3:4];
                    obj.x = obj.State(obj.state_to_x);
                    obj.sigma_r = sigma_r;
                    obj.R = sigma_r * [1 0; 0 1];
                    obj.sigma_q = sigma_q;
                    obj.Q = sigma_q * eye(obj.dims);
                    %obj.Q(3,3) = sigma_q;
                    %obj.Q(6,6) = sigma_q;
                    obj.m = m;
                    obj.time = time;
                    % States = [x x_vel f_x y y_vel f_y];
                    obj.A = [0 1 0 0 0 0; 0 0 1/m 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1/m; 0 0 0 0 0 1];
                    obj.H = [1 0 0 0 0 0; 0 0 0 1 0 0];
                    obj.A_discr = time*obj.A + eye(obj.dims);
                case 'speed'
                    obj.dims = 4;
                    obj.type = type;
                    obj.P = sigma_p*eye(obj.dims);
                    obj.P(1,1) = 0;
                    obj.P(3,3) = 0;
                    obj.sigma_p = sigma_p;
                    if(length(init_x)==2)
                        obj.State = [init_x(1);0;init_x(2);0];
                    else
                        obj.State = [init_x(1);init_x(3);init_x(2);init_x(4)];
                    end
                    obj.state_to_x = [1:2:3];
                    obj.x = obj.State(obj.state_to_x);
                    obj.sigma_r = sigma_r;
                    obj.R = sigma_r * [1 0; 0 1];
                    obj.sigma_q = sigma_q;
                    obj.Q = eye(obj.dims);
                    obj.Q(1,1) = sigma_q(1);
                    obj.Q(3,3) = sigma_q(1);
                    obj.Q(2,2) = sigma_q(2);
                    obj.Q(4,4) = sigma_q(2);
                    obj.m = m;
                    obj.time = time;
                    % States = [x x_vel y y_vel];
                    obj.A = [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];
                    obj.H = [1 0 0 0; 0 0 1 0];
                    obj.A_discr = time*obj.A + eye(obj.dims);
                case 'acceleration'
                    obj.dims = 6;
                    obj.type = type;
                    obj.P = sigma_p*eye(obj.dims);
                    obj.P(1,1) = 0;
                    obj.P(4,4) = 0;
                    obj.sigma_p = sigma_p;
                    obj.State = [init_x(1);0;0;init_x(2);0;0];
                    obj.state_to_x = [1:3:4];
                    obj.x = obj.State(obj.state_to_x);
                    obj.sigma_r = sigma_r;
                    obj.R = sigma_r * [1 0; 0 1];
                    obj.sigma_q = sigma_q;
                    obj.Q = eye(obj.dims);
                    obj.Q(1,1) = sigma_q(1);
                    obj.Q(4,4) = sigma_q(1);
                    obj.Q(2,2) = sigma_q(2);
                    obj.Q(5,5) = sigma_q(2);
                    obj.Q(3,3) = sigma_q(3);
                    obj.Q(6,6) = sigma_q(3);
                    obj.m = m;
                    obj.time = time;
                    % States = [x x_vel x_acc y y_vel y_acc];
                    obj.A = [0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 0 0 0 0 0 1];
                    obj.H = [1 0 0 0 0 0; 0 0 0 1 0 0];
                    obj.A_discr = time*obj.A + eye(obj.dims);
                case 'custom'
                otherwise
                    disp(['Wrong initialization with non-existent type ' type ' try force, speed, acceleration or custom instead.']);
            end
        end
        
        function [obj,x,P] = Estimate(obj)
            obj.State = obj.A_discr * obj.State;
            obj.P = obj.A_discr*obj.P*obj.A_discr' + obj.Q;
            obj.x = obj.State(obj.state_to_x); 
            x = obj.x;
            P = obj.P;
        end
        
        function obj = Update(obj,x_meas,m)
            obj.m = m;
            if (strcmp(obj.type,'force'))
                obj.A = [0 1 0 0 0 0; 0 0 1/m 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1/m; 0 0 0 0 0 1];
                obj.A_discr = obj.time*obj.A + eye(obj.dims);
            end
            K = obj.P * obj.H' * inv(obj.H * obj.P * obj.H' + obj.R);
            
            State = obj.State/obj.m;
            x_meas_v2 = x_meas/obj.m;
            %obj.State = obj.State + K * (x_meas - obj.H * obj.State);
            State = State + K * (x_meas_v2 - obj.H * State);
            obj.State = State*obj.m;
            obj.P = (eye(obj.dims) - K * obj.H) * obj.P;
            obj.x = obj.State(obj.state_to_x);
        end
    end
end

