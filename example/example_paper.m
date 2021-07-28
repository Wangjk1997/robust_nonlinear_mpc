clear all
addpath('../src/')
addpath('../src/utils/')
mysys = nonlinearDiscreteSystem(2,1);
x = ones(2,1);
xHistory = [];
n_horizon = 10;
Xc = {[],[]};
Uc = {[-0.475,0.475]};
nominalMPC = nonlinearMPC(mysys, Xc, Uc, n_horizon, eye(2,2), 0.2);
nominalMPC = nominalMPC.add_initial_constraint(x);
[x,u] = nominalMPC.solve()
