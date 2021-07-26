classdef nonlinearMPC
    properties
        sys; %system
        Xc; Uc; % constraints set for statespace and input space
        Q; R;
        n_horizon; % prediction horizon
        n_opt; % dim. of optimization parameter s :=[x(0), .., x(n), u(0)....u(n-1)]
        H; % positive definite matrix for objective function V(s) = s'*H*s
        Cineq_A;
        Cineq_b;
        Cinit_A;
        Cinit_b;
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
                state_in = s(obj.sys.nx * (i - 1) + 1:obj.sys.nx * i);
                state_out = s(obj.sys.nx * i + 1:obj.sys.nx * (i + 1));
                input_start_index = obj.sys.nx * (obj.n_horizon + 1);
                input = s(input_start_index + obj.sys.nu * (i - 1) + 1: input_start_index + obj.sys.nu * i);
                ceq(obj.sys.nx * (i - 1) + 1:obj.sys.nx * i) = obj.sys.propagate(state_in, input) - state_out;
            end
        end
        
        function [C_ineq_A, C_ineq_b] = construct_ineq_constraint(obj, Xc, Uc)
           number_constraint_Xc = 0;
           index_Xc = [];
           for i = 1:size(Xc,1)
               if(~isempty(Xc{i}))
                   number_constraint_Xc = number_constraint_Xc + 1;
                   index_Xc = [index_Xc, i];
               end
           end
           number_constraint_Uc = 0;
           index_Uc = [];
           for i = 1:size(Uc,1)
               if(~isempty(Uc{i}))
                   number_constraint_Uc = number_constraint_Uc + 1;
                   index_Uc = [index_Uc, i];
               end
           end
           C_ineq_A_Xc = zeros(2 * (number_constraint_Xc * (obj.n_horizon + 1)),obj.n_opt);
           C_ineq_A_Uc = zeros(2 * number_constraint_Uc * obj.n_horizon, obj.n_opt);
           C_ineq_b_Xc = zeros(2 * (number_constraint_Xc * (obj.n_horizon + 1)), 1);
           C_ineq_b_Uc = zeros(2 * number_constraint_Uc * obj.n_horizon, 1);
           for i = 1:obj.n_horizon + 1
               start_index_row = 2 * number_constraint_Xc * (i-1);
               start_index_col = obj.sys.nx * (i-1);
               for j = 1:number_constraint_Xc
                   boundary = Xc{index_Xc(j)};
                   lb = boundary(1);
                   ub = boundary(2);
                   if lb > 0
                       tmp_lb_A = -1/lb;
                       tmp_lb_b = -1;
                   elseif lb < 0
                       tmp_lb_A = 1/lb;
                       tmp_lb_b = 1;
                   else
                       tmp_lb_A = -1;
                       tmp_lb_b = 0;
                   end
                   if ub > 0
                       tmp_ub_A = 1/ub;
                       tmp_ub_b = 1;
                   elseif ub < 0
                       tmp_ub_A = -1/ub;
                       tmp_ub_b = -1;
                   else
                       tmp_ub_A = 1;
                       tmp_ub_b = 0;
                   end
                   C_ineq_A_Xc(start_index_row + 2*(j-1) + 1, start_index_col + index_Xc(j)) = tmp_lb_A;
                   C_ineq_b_Xc(start_index_row + 2*(j-1) + 1, 1) = tmp_lb_b;
                   C_ineq_A_Xc(start_index_row + 2*(j-1) + 2, start_index_col + index_Xc(j)) = tmp_ub_A;
                   C_ineq_b_Xc(start_index_row + 2*(j-1) + 2, 1) = tmp_ub_b;
               end
           end
           for i = 1:obj.n_horizon
               start_index_row = 2 * number_constraint_Uc * (i-1);
               start_index_col = obj.sys.nx * (obj.n_horizon + 1) + obj.sys.nu * (i-1);
               for j = 1:number_constraint_Uc
                   boundary = Uc{index_Uc(j)};
                   lb = boundary(1);
                   ub = boundary(2);
                   if lb > 0
                       tmp_lb_A = -1/lb;
                       tmp_lb_b = -1;
                   elseif lb < 0
                       tmp_lb_A = 1/lb;
                       tmp_lb_b = 1;
                   else
                       tmp_lb_A = -1;
                       tmp_lb_b = 0;
                   end
                   if ub > 0
                       tmp_ub_A = 1/ub;
                       tmp_ub_b = 1;
                   elseif ub < 0
                       tmp_ub_A = -1/ub;
                       tmp_ub_b = -1;
                   else
                       tmp_ub_A = 1;
                       tmp_ub_b = 0;
                   end
                   C_ineq_A_Uc(start_index_row + 2*(j-1) + 1, start_index_col + index_Uc(j)) = tmp_lb_A;
                   C_ineq_b_Uc(start_index_row + 2*(j-1) + 1, 1) = tmp_lb_b;
                   C_ineq_A_Uc(start_index_row + 2*(j-1) + 2, start_index_col + index_Uc(j)) = tmp_ub_A;
                   C_ineq_b_Uc(start_index_row + 2*(j-1) + 2, 1) = tmp_ub_b;
               end
           end
           C_ineq_A = [C_ineq_A_Xc; C_ineq_A_Uc];
           C_ineq_b = [C_ineq_b_Xc; C_ineq_b_Uc];
        end
        
        function obj = add_initial_constraint(obj, x_init)
            idx_start = 1;
            idx_end = obj.sys.nx;
            obj.Cinit_A = zeros(obj.sys.nx, obj.n_opt);
            obj.Cinit_A(:, idx_start:idx_end) = eye(obj.sys.nx);
            obj.Cinit_b = x_init;
        end
        
        
        function [x_seq, u_seq] = solve(obj)
            options = optimoptions('fmincon', 'Display', 'none');
            cost_function = @(s)obj.costfunction(s);
            A = obj.Cineq_A;
            b = obj.Cineq_b;
            Aeq = obj.Cinit_A;
            beq = obj.Cinit_b;
            nonlcon = @(s)obj.dynamic_nonlinear_constraint(s);
            % assume u = 0 is the initial searching point.
            s0 = zeros(obj.n_opt,1);
            s0(1:obj.sys.nx) = obj.Cinit_b;
            x = obj.Cinit_b;
            for i = 1:obj.n_horizon
                u_next = zeros(obj.sys.nu,1);
                x = mysys.propagate(x, u_next);
                s0(obj.sys.nx*(i)+1:obj.sys.nx*(i)+1) = x;
            end
            [var_optim] = fmincon(cost_function, s0, A, b, Aeq, beq,[],[],nonlcon, options);
            x_seq = reshape(var_optim(1:obj.sys.nx*(obj.n_horizon+1)), obj.sys.nx, obj.n_horizon+1);
            u_seq = reshape(var_optim(obj.sys.nx*(obj.n_horizon+1)+1:obj.n_opt), obj.sys.nu, obj.n_horizon);
        end
    end
    
end