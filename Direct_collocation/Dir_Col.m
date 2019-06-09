% Code: Direct Collocation to solve the optimal control problem

clc
clear 
close all 

% load an optimised solution to initialise
load('DirColInitial4')

% mesh size
gridN = 300;
tic

% objective function
obj_fun = @cc_obj_fun;

% The initial parameter guess vector = [t_f; x1; x2; x3; x4 ; u]
% A ramp
% x0 = [120; ones(gridN,1); linspace(0.2,0,gridN)'; linspace(-1.5,0,gridN)';...
%     linspace(1,0,gridN)'; linspace(3,0,gridN)'];

% load a partially optimised solution to initialise
load('DirColInitial4')
x0 = [120; x1I'; x2I'; x3I'; x4I' ; x5I'];

% No linear inequality or equality constraints
A = [];
b = [];
Aeq = [];
Beq = [];

% Lower bound and upper bound
lb = [0;    ones(gridN , 1) * (0 - 1);  zeros(gridN, 1); ones(gridN , 1) * -1.5;
    zeros(gridN, 1); zeros(gridN, 1)];
ub = [Inf;  ones(gridN , 1) * inf;   ones(gridN, 1) *inf ; ones(gridN , 1) * inf;
    ones(gridN , 1)*inf ; ones(gridN, 1)*inf];

% Options for fmincon
options = optimoptions(@fmincon, 'TolFun', 0.00000001, 'MaxIter', 10000, ...
                       'MaxFunEvals', 1000000, 'Display', 'iter', ...
                       'DiffMinChange', 0.001, 'Algorithm', 'sqp');
                   
% Solve for the best simulation time + control input
[optimal,fval] = fmincon(obj_fun, x0, A, b, Aeq, Beq, lb, ub, ...
              @cc_cons, options);

%% Plotting the solution         
% Discretize the times
sim_time = optimal(1);
delta_time = sim_time / gridN;
times = 0 : delta_time : sim_time - delta_time;

% Get the states and control out of the vector
x1      = optimal(2             : 1 + gridN);
x2      = optimal(2 + gridN   : 1 + gridN * 2);
x3      = optimal(2 + gridN *2  : 1 + gridN * 3);
x4      = optimal(2 + gridN *3  : 1 + gridN * 4);
u       = optimal(2 + gridN *4  : end);

%shift back from the origin
x1 = x1 + 1;
x3 = x3 + 1.65;

% Make the plots
figure();
plot(times, x1,'o',times, x2,'o',times, x3,'o',times, x4,'o',times,u,'o');
legend('Normal Cell','Tumor Cell','Immune Cell','Drug Concentration','Control Input')
title('Optimal Output');
xlabel('Time (s)');
ylabel('Cell Numbers  (m/s^2)');
grid on

disp(sprintf('Finished in %f seconds', toc));


%% objective function
function [cost,grad] = cc_obj_fun(X_opt)
% X_opt input is in the form of [t_final;x1;x2;x3;x4;u]

%number of grids used in discretisation
gridN = (length(X_opt) - 1)/5;

% time discretisation
t_f = X_opt(1);
t = linspace(0,1,gridN)'*t_f;
dt = t(2) - t(1);

% number of states
n_s = 4;

% reshape the states into [x1;x2;x3;x4]
x = X_opt(2:length(X_opt) - gridN);
x = reshape(x,gridN,n_s);

% get control
u = X_opt(length(X_opt) - gridN + 1:end);

% calculate cost using trapzoidal integration
Q =[0,100,0,0.1]';
y = 1/2*(x.^2*Q + u.^2);
cost = trapz(t,y);

% get gradient 
grad_states = [x(1,:)*diag(Q)*dt; x(2:end-1,:)*diag(Q)*2*dt; x(end,:)*diag(Q)*dt];
grad_control = [u(1)*dt; u(2:end-1)*2*dt; u(end)*dt];

grad = ones(size(X_opt));
grad(1) = 0;
grad(2:end) = [grad_states(:);grad_control];

end

%% constraints of direct collocation
function [C , Ceq] = cc_cons(X_opt)


%number of grids used in discretisation
gridN = (length(X_opt) - 1)/5;

% time discretisation
t_f = X_opt(1);
t = linspace(0,1,gridN)'*t_f;
dt = t(2) - t(1);

% number of states
n_s = 4;

% reshape the states into [x1;x2;x3;x4]
x = X_opt(2:length(X_opt) - gridN);
x = reshape(x,gridN,n_s);

x_ll = x(1:end-1,:);
x_rr = x(2:end,:);

% get control
u = X_opt(length(X_opt) - gridN + 1:end);

% get f(x,u)_k and f(x,u)_k+1
x_dot = cc_fun(t,x,u);
x_dot_ll = x_dot(1:end-1,:);
x_dot_rr = x_dot(2:end,:);

% get collocation points time control and states
t_c = (t(1:end-1) + t(2:end))/2;
u_c = (u(1:end-1) + u(2:end))/2;
x_c = 1/2 * (x_ll + x_rr) + dt/8*(x_dot_ll - x_dot_rr);
x_dot_c = cc_fun(t_c,x_c,u_c);

% the constraint of the dynamical system 
Ceq = x_ll - x_rr + dt/6*(x_dot_ll + x_dot_rr + 4*x_dot_c);
Ceq = Ceq(:);

% initial and final constraints
x1_0 = 0;
x2_0 = 0.2;
x3_0 = -1.50;
x4_0 = 0;


% input into equality constraints
Ceq = [Ceq; x(1,1) - x1_0; x(1,2) - x2_0; x(1,3) - x3_0; x(1,4) - x4_0;
    x(end,1); x(end,2); x(end,3); x(end,4)];
C = [];

end

%% calculate the rates
function [x_dot] = cc_fun(t,x,u)

x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);
x4 = x(:,4);

% parameters
a1 = 0.2;
a2 = 0.3;
a3 = 0.1;
b1 = 1;
b2 = 1;
alpha = 0.3;
c1 = 1;
c2 = 0.5;
c3 = 1;
c4 = 1;
d1 = 0.2;
d2 = 1;
r1 = 1.5;
r2 = 1;
s = 0.33;
rho = 0.01;


% non-linear dynamics shifted
x1_dot = -r2*x1.*(1+b2*x1) - c4*x2.*x1 - c4/b2*x2 -a3*x4.*x1 - a3/b2*x4;
x2_dot = r1*x2.*(1-b1*x2) - (c2*s/d1 + c3/b2)*x2 - c3*x2.*x1 - c2*x2.*x3 - a2*x4.*x2;
x3_dot =  - c1*s/d1*x2 - d1*x3 - a1*s/d1*x4 + rho*s/d1*x2./(alpha + x2) + rho*x2.*x3./(alpha + x2) ...
    - c1*x3.*x2 - a1*x4.*x3;
x4_dot = u - d2*x4;


x_dot = [x1_dot,x2_dot,x3_dot,x4_dot];

end
