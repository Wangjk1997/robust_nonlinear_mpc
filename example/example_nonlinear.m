clear all
addpath('../src/')
addpath('../src/utils/')
mysys = nonlinearSystem(6,2,1/120,1/1000);
x = zeros(6,1);
s = [x];
xHistory = [];
n_horizon = 5;
for i = 1:n_horizon
    u_next = [0.01;0];
    x = mysys.propagate(x, u_next);
    xHistory = [xHistory,x];
    s = [s;x];
end
for i = 1:n_horizon
    s = [s;u_next];
end
plot(xHistory(1,:))

Xc = {[];[0, 0.3];[];[-0.3,0.3];[-pi/3,pi/3];[]};
Uc= {[-0.02, 0.02];[-0.02, 0.02]};


nominalMPC = nonlinearMPC(mysys, Xc, Uc, n_horizon, eye(6,6), 0.5*eye(2,2))
H = nominalMPC.H;
[c1,c2] = nominalMPC.dynamic_nonlinear_constraint(s);