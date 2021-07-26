clear all
addpath('../src/')
addpath('../src/utils/')
mysys = nonlinearSystem(6,2,1/120,1/1000);
x = zeros(6,1);
xHistory = [];
for i = 1:(10*120)
    x = mysys.propagate(x,[0.01;0]);
    xHistory = [xHistory,x];
end
plot(xHistory(5,:))
% 
% nominalMPC = nonlinearMPC(mysys, 1, 2, 2, eye(6,6), 0.5*eye(2,2))
% H = nominalMPC.H;
% [c1,c2] = nominalMPC.dynamic_nonlinear_constraint(s)