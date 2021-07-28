classdef nonlinearDiscreteSystem
    properties
        nx; nu; % dim of state space and input space
    end
    methods
        function obj = nonlinearDiscreteSystem(nx, nu)
            obj.nx = nx;
            obj.nu = nu;
        end
        
        function x_new = propagate(obj, x, u) 
            % put your nonlinear system state function here
            % initialize non-inputs
            x1 = x(1);
            x2 = x(2);
            x_new = zeros(2,1);
            x_new(1) = x2;
            x_new(2) = sin(x1) + u;
        end
    end
end