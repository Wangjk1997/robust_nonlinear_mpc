example = exampleClass(eye(2,2));
example.fun([1;1]);
fun = @(x)example.fun(x);
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = @(x)example.constraint(x);
x0 = [0;0];
options = optimoptions('fmincon','Display', 'none');
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);