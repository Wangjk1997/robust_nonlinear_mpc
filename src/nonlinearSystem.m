classdef nonlinearSystem
    properties
        nx; nu; % dim of state space and input space
        Dt; %Samlping time
        dt; %Step time
    end
    methods
        function obj = nonlinearSystem(nx, nu, Dt, dt)
            obj.nx = nx;
            obj.nu = nu;
            obj.Dt = Dt;
            obj.dt = dt;
        end
        
        function x_new = propagate(obj, x, u) 
            % put your nonlinear system state function here
            % initialize non-inputs
            mx = 0.1708;
            mz = 0.1794;
            DvxCB = 0.0125;
            DvzCB = 0.0480;
            rz = 0.09705; %r_(z,g/b)^b
            mRB = 0.1249;
            Icgy = 0.005821;
            DomegaCB = 0.000862;
            g = 9.8;
            mAx = mx - mRB;
            mAz = mz - mRB;
            T = 0.26; %r_(z,t/b)^b
            
            for i = 1:(obj.Dt/obj.dt)
                % state vector
                x1 = x(1); % px_g/n_b
                x2 = x(2); % vx_g/n_b
                x3 = x(3); % pz_g/n_b
                x4 = x(4); % vx_g/n_b
                x5 = x(5); % theta
                x6 = x(6); % omega
            
                % input vector
                u1 = u(1);
                u2 = u(2);
                u3 = (T - rz) * u(1);
            
                % state function
                dxdt = zeros(6,1);
                dxdt(1) = cos(x5) * x2 + sin(x5) * x4;
                dxdt(2) = -(DvxCB*Icgy*x2 - Icgy*u1 - mAx*rz*u3 + Icgy*mAz*x4*x6 + Icgy*mRB*x4*x6 - DvxCB*mAx*rz^2*x2 + DvxCB*mAx*rz^3*x6 + mAx^2*rz*x2*x4 - DvxCB*Icgy*rz*x6 - mAx^2*rz^2*x4*x6 + DomegaCB*mAx*rz*x6 + g*mAx*mRB*rz^2*sin(x5) - mAx*mAz*rz*x2*x4)/(- mAx^2*rz^2 + Icgy*mAx + Icgy*mRB);
                dxdt(3) = -sin(x5) * x2 + cos(x5) * x4;
                dxdt(4) = (u2 - DvzCB*x4 - mAx*rz*x6^2 + mAx*x2*x6 + mRB*x2*x6)/(mAz + mRB);
                dxdt(5) = x6;
                dxdt(6) = -(mAx^2*x2*x4 - mRB*u3 - mAx*u3 + DomegaCB*mAx*x6 + DomegaCB*mRB*x6 - mAx*rz*u1 - mAx*mAz*x2*x4 + mAx*mRB*x2*x4 - mAz*mRB*x2*x4 + DvxCB*mRB*rz^2*x6 - mAx^2*rz*x4*x6 + g*mRB^2*rz*sin(x5) - DvxCB*mRB*rz*x2 + mAx*mAz*rz*x4*x6 + g*mAx*mRB*rz*sin(x5))/(- mAx^2*rz^2 + Icgy*mAx + Icgy*mRB);
                
                x = x + obj.dt * dxdt;
            end
            
            % get new state vector
            x_new = x;
        end
    end
end