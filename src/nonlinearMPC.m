classdef nonlinearMPC
    properties
        sys; %system
        Xc; Uc; % constraints set for statespace and input space
        Q; R;
        x_min; x_max; % lower and upper bound of Xc
        n_horizon; % prediction horizon
        n_opt; % dim. of optimization parameter s :=[x(0), .., x(n), u(0)....u(n-1)]
        H; % positive definite matrix for objective function V(s) = s'*H*s
        Cineq_A;
        Cineq_b;
    end
    methods(Access = public)
        function obj = nonlinearMPC(sys, Xc, Uc, N, Q, R)
            obj.sys = sys;
            obj.Xc = Xc;
            obj.Uc = Uc;
            obj.n_horizon = N;
            obj.n_opt = obj.sys.nx*(obj.n_horizon+1)+obj.sys.nu*obj.n_horizon;
            obj.Q = Q;
            obj.R = R;
            obj.H = obj.construct_costfunction_matrix();
            [obj.Cineq_A, obj.Cineq_b] = obj.construct_ineq_constraint(Xc,Uc);
        end
    end
    methods(Access = public)
        function H = construct_costfunction_matrix(obj)
            Q_block = [];
            R_block = [];
            for itr=1:obj.n_horizon
                Q_block = blkdiag(Q_block, obj.Q);
                R_block = blkdiag(R_block, obj.R);
            end
            H = blkdiag(Q_block, zeros(obj.sys.nx, obj.sys.nx), R_block); % here we assume V_f(x)=0;
        end
        
        function value = costfunction(obj, s)
            value = s'* obj.H * s;
        end
        
        function [c, ceq] = dynamic_nonlinear_constraint(obj, s)
            c = [];
            ceq = zeros(obj.sys.nx * obj.n_horizon, 1);
            for i = 1:obj.n_horizon
                state_in = s(obj.sys.nx * (i - 1) + 1:obj.sys.nx * i)
                state_out = s(obj.sys.nx * i + 1:obj.sys.nx * (i + 1))
                input_start_index = obj.sys.nx * (obj.n_horizon + 1);
                input = s(input_start_index + obj.sys.nu * (i - 1) + 1: input_start_index + obj.sys.nu * i)
                ceq(obj.sys.nx * (i - 1) + 1:obj.sys.nx * i) = obj.sys.propagate(state_in, input) - state_out;
            end
        end
        
        function [C_ineq_A, C_ineq_b] = construct_ineq_constraint(obj, Xc, Uc)
            % compute C_ineq
            [F, G, nc] = convert_Poly2Mat(Xc, Uc);
            F_block = [];
            G_block = [];
            for itr = 1:obj.n_horizon
                G_block = blkdiag(G_block, G);
            end
            for itr = 1:obj.n_horizon+1
                F_block = blkdiag(F_block, F);
            end
            C_ineq_A = [F_block, [G_block; zeros(nc, obj.sys.nu*obj.n_horizon)]];
            nc_total = size(C_ineq_A, 1);
            C_ineq_b = ones(nc_total, 1);
        end
%         function [x_seq, u_seq] = solve(obj)
%             
%         end
    end
    
end