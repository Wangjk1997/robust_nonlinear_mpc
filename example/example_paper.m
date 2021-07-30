clear all
addpath('../src/')
addpath('../src/utils/')
mysys = nonlinearDiscreteSystem(2,1);
x = [0.8;1];
xHistory = x;
uHistory = [];
n_horizon = 6;
Xc = {[],[]};
Uc = {[-0.475,0.475]};
nominalMPC = nonlinearMPC(mysys, Xc, Uc, n_horizon, eye(2,2), 0.5);
nominalMPC = nominalMPC.set_reference([0.5;0.5]);
for i = 1:50
    nominalMPC = nominalMPC.add_initial_constraint(x);
    [x_seq,u_seq] = nominalMPC.solve();
    x = mysys.propagate(x, u_seq(1));
    xHistory = [xHistory,x];
    uHistory = [uHistory,u_seq(1)];
end
figure(1)
plot(0:10,xHistory(1,1:11))
xlabel("Time")
ylabel("z_1")
x0=100; y0=100;
width=400; height=220;
set(gcf,'units','points','position',[x0,y0,width,height])
% print('-dpng','-r300','x1.png')
figure(2)
plot(0:10,xHistory(2,1:11))
xlabel("Time")
ylabel("z_2")
set(gcf,'units','points','position',[x0,y0,width,height])
% print('-dpng','-r300','x2.png')
figure(3)
plot(0:9,uHistory(1:10))
xlim([0,10])
xlabel("Time")
ylabel("u")
set(gcf,'units','points','position',[x0,y0,width,height])
% print('-dpng','-r300','u.png')
figure(4)
xlabel("Time")
ylabel("z")
hold on
p1 = plot(0:10,xHistory(1,1:11));
p2 = plot(0:10,xHistory(2,1:11));
legend([p1 p2],{'z_1', 'z_2'})
set(gcf,'units','points','position',[x0,y0,width,height])
% print('-dpng','-r300','z.png')
figure(5)
plot(xHistory(1,1:11),xHistory(2,1:11))
xlabel("z_1")
ylabel("z_2")
axis equal
