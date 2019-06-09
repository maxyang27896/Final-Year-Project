% Code: Optimising the dynamic control variables using Genetic Algorithm

clear
clc
close all

% Define ploblem
ObjectiveFunction = @cancer_opt_fun;
nvars = 10;    % Number of variables
A = []; 
B = [];
Aeq = [];
Beq = [];

% Hard variable constraints
LB = [-inf,-inf,-inf,  0,  0,   0, 0,    0, 0,    0];   % Lower bound
UB = [ inf, inf, inf,inf,  inf, inf, 1,  inf, 1,  inf];  % Upper bound

% Initial population constrains 
initialLB = [-200,-200,-200,   0, 0, 0,    0,   0,    0,    0];
initialUB = [ 200, 200, 200, 200, 1, 0.1,  0.5,   5,  0.5,    5];
nonlcon = [];
intcon = []; %set all var to be integers 
options = optimoptions('ga','ConstraintTolerance',1e-6','PlotFcn',{@gaplotbestf,@gaplotdistance},...
        'InitialPopulationRange',[initialLB;initialUB],'FitnessScalingFcn',...
        @fitscalingprop,'SelectionFcn',@selectionremainder,...
        'EliteCount',20,'CrossoverFraction',0.9,'CrossoverFcn',@crossoverheuristic,...
        'UseParallel',true);
    
% solve 
[x,fval,exitflag,output,populations,scores] = ga(ObjectiveFunction,nvars,...
    A,B,Aeq,Beq,LB,UB,nonlcon,intcon,options);


%% run the dynamic control law results and plot the solution
a1 = 0.2;
a2 = 0.3;
a3 = 0.1;
b1 = 1;
b2 = 1;
c1 = 1; 
c2 = 0.5;
c3 = 1;
c4 = 1;
d1 = 0.2;
d2 = 1;
r1 = 1.5;
r2 = 1;
s = 0.33;
alpha = 0.3;
rho = 0.01;

tspan = [0,350];

% initial conditions
N0 = 1;
T0 = 0.2;
I0 = 0.15;
M0 = 0;    
J0 = 0;    
C0 = 0;     

% origin shift
x1_0 = N0 - 1/b2;
x2_0 = T0;
x3_0 = I0 - s/d1;
x4_0 = M0;


% paramters taken from GA output
xi1_0 = x(1);
xi2_0 = x(2);
xi3_0 = x(3);
xi4_0 = x(4);
k = x(5);
r_constant = x(6);
gamma1 = x(7);
gamma2 = x(8);
gamma3 = x(9);
gamma4 = x(10);

% input paramters and initial conditions to simulate system
y0_DCL = [x1_0;x2_0;x3_0;x4_0;xi1_0;xi2_0;xi3_0;xi4_0;J0;C0];
opt    = odeset('Events', @myEvent);
[t2,y2] = ode45(@(t,y) nonlinear_rhs_opt(t,y,k,r_constant,gamma1,gamma2,gamma3,...
    gamma4), tspan, y0_DCL, opt);

% getting the controls and additional cost
u_dynamic = zeros(1,length(y2));
c_x_xi = zeros(1,length(y2));

for i = 1:length(y2)
    [~,u_dynamic(:,i),c_x_xi(:,i)] = nonlinear_rhs_opt(t2(i),y2(i,:),k,r_constant,...
        gamma1,gamma2,gamma3,gamma4);
end
y2(end,9:10)

% shift back the origin
y2(:,1) = y2(:,1) + 1/b2;
y2(:,3) = y2(:,3) + s/d1;

% plot states
figure
plot(t2,y2(:,1:4),'-',t2,u_dynamic,'-','linewidth',3)
hold on
plot(tspan,[0.75,0.75],'--k','linewidth',2)
legend('Normal Cell','Tumor Cell','Immune Cell','Drug Dose','Control Input','Normal Cell Limit')
grid on
xlabel('Time (days)')
ylabel('Normalised Fraction of Cell')

% plot the augmented states 
figure
plot(t2,y2(:,5:8),'-','linewidth',3)
legend('\xi_1','\xi_2','\xi_3','\xi_4')
grid on
xlabel('Time (days)')
ylabel('Augmented States (\xi)')

% plot the cost
figure
plot(t2,y2(:,9:10),'-','linewidth',3)
legend('Tumor Cost','Approximated Cost')
grid on
xlabel('Time (days)')
ylabel('Normalised Fraction of Cell')


%% main code end---------------------------------------------------------------------------------------------------------
% Objective Function for ga
function cost = cancer_opt_fun(x)

% cancer model constants 
a1 = 0.2;
a2 = 0.3;
a3 = 0.1;
b1 = 1;
b2 = 1;
c1 = 1;
c2 = 0.5;
c3 = 1;
c4 = 1;
d1 = 0.2;
d2 = 1;
r1 = 1.5;
r2 = 1;
s = 0.33;
alpha = 0.3;
rho = 0.01;

%% simulate trajectory: 
%state = [Normal Cell, Tumor Cell, Immune Cell, Drug, Cost, Additional Cost]

tspan = [0,50];

% initial condtions
N0 = 1;
T0 = 0.2;
I0 = 0.15;
M0 = 0;     
J0 = 0;     
C0 = 0;     

% origin shift
x1_0 = N0 - 1/b2;
x2_0 = T0;
x3_0 = I0 - s/d1;
x4_0 = M0;

% augmented states 
xi1_0 = x(1);
xi2_0 = x(2);
xi3_0 = x(3);
xi4_0 = x(4);

% dynamic control paramters 
k = x(5);
r_constant  = x(6);

%% run the nonlinear dynamic control law
y0_DCL = [x1_0;x2_0;x3_0;x4_0;xi1_0;xi2_0;xi3_0;xi4_0;J0;C0];

% algebraic P solution tuning parameters
gamma_1 = x(7);
gamma_2 = x(8);
gamma_3 = x(9);
gamma_4 = x(10);

opt    = odeset('Events', @myEvent);
[t,y] = ode45(@(t,y) nonlinear_rhs_opt(t,y,k,r_constant,gamma_1,gamma_2,gamma_3...
    ,gamma_4), tspan, y0_DCL, opt);

%% define cost 
cost = y(end,10);

% if inequality not satisfied, set arbitrary large cost
if isnan(cost) || isinf(cost)
    cost = 10000;
    return
end

end


%% rates function
function [y_dot,u,c_x_xi] = nonlinear_rhs_opt(t,y,k,r_constant,gamma1,gamma2_bar,gamma3,gamma4_bar)
% define ODE dynamics here to be used in solver, states x = [N,I,T,M]' 

% input y-vector 
x1 = y(1);
x2 = y(2);
x3 = y(3); 
x4 = y(4);

x = [x1;x2;x3;x4];

% augmented states
xi1 = y(5);
xi2 = y(6);
xi3 = y(7);
xi4 = y(8);

xi = [xi1;xi2;xi3;xi4];

% cancer model constants
a1 = 0.2;
a2 = 0.3;
a3 = 0.1;
b1 = 1;
b2 = 1;
c1 = 1;
c2 = 0.5;
c3 = 1;
c4 = 1;
d1 = 0.2;
d2 = 1;
r1 = 1.5;
r2 = 1;
s = 0.33;
alpha = 0.3;
rho = 0.01;

% penalty parameters
Q11 = 0;
Q22 = 100;
Q33 = 0;
Q44 = 0.1;
Q = diag([Q11,Q22,Q33,Q44]);

% dynamic control law parameter
R = eye(4)*r_constant;


%%  the linear parametrisation F(x) and control matrix B(x)
F11 = -r2*(x1*b2+1) - c4*x2 - a3*x4;
F12 = -c4/b2;
F13 = 0;
F14 = -a3/b2;

F21 = 0;
F22 = r1*(1-b1*x2) - (c2*s/d1 + c3/b2) - c3*x1 - c2*x3 - a2*x4;
F23 = 0;
F24 = 0;

F31 = 0;
F32 = -c1*s/d1 + rho*s/d1/(alpha + x2);
F33 = -d1 + rho*x2/(alpha+x2) - c1*x2 - a1*x4;
F34 = -a1*s/d1;

F41 = 0;
F42 = 0;
F43 = 0;
F44 = -d2;

F = [F11,F12,F13,F14;F21,F22,F23,F24;F31,F32,F33,F34;F41,F42,F43,F44];

B = [0;0;0;1];


%% algebraic P solutions paramter derivation
gamma4 = -sqrt(Q44 + 1) + sqrt(Q44+1+2.7225*Q33) + gamma4_bar;
K4 = -1 + sqrt(Q44+1) + gamma4;
D4 = K4^2 + 2*K4 - Q44;

alpha3 = 250000/1089*D4;

% average between 0 and 1 
K3 =  alpha3*(gamma3*(1/25 - 1/25*sqrt(1-Q33/alpha3*25^2)) +...
    (1-gamma3)*(1/25 + 1/25*sqrt(1-Q33/alpha3*25^2)));
D3 = 2/25*K3 - Q33 - 1089/250000*K3^2*1/D4;


gamma2 = (1+0.0021*K3^2/(D4*D3))^2/(1/Q11 - 1/(100*D4)) + gamma2_bar;
K2 = 800/169*(Q22 + 0.1018/(D3)*K3^2 + gamma2);
D2 = 169/800*K2 - Q22 - 0.1018*K3^2/D3;

alpha1 = 1/(1/(100*D4) + (1+0.0021*K3^2/(D4*D3))^2/D2);

% an average between 0 and 1 
K1 = alpha1*(gamma1*(1-sqrt(1-Q11/alpha1)) + (1-gamma1)*(1+sqrt(1-Q11/alpha1)));
D1 = -K1^2*(1/alpha1) + 2*K1 - Q11;

%% algebraic P solution  
P_xi_jacobian_x = [ K1*x1,             K1*x1,         0,   (K1*x1)/10;
 K2*x2,                                (3*K2*x2)/2, (K2*x2)/2, (3*K2*x2)/10;
     0, K3*x3*(xi2/(100*(xi2 + 3/10)^2) - 1/(100*(xi2 + 3/10)) + 1),         0,    (K3*x3)/5;
     0,                                 0,         0,            0];

P44_xi = K4;
P33_xi = K3*(xi2 + xi4/5 - xi2/(100*(xi2 + 3/10)) + 1/5) ;
P22_xi = K2*(xi1 + (3*xi2)/2 + xi3/2 + (3*xi4)/10 + 13/40);
P11_xi = K1*(xi1 + xi2 + xi4/10 + 1);

P_xi = diag([P11_xi,P22_xi,P33_xi,P44_xi]);


%% value function derivatives 
dv_dx = x'* P_xi +  (x - xi)'*R;
dv_dxi =  1/2 * x'* P_xi_jacobian_x - ( x - xi)'*R;


%% ==========================================================================================
% input dynamic control law
u = - (B' * dv_dx');

%% check the inequlaity conditions and stop ode45 if not satisfied
c_x_xi = dv_dx*F*x + 1/2*x' * Q * x - ...
    1/2*dv_dx*B*B'*dv_dx'  - k* dv_dxi * dv_dxi';

if c_x_xi > 0 
    y_dot = inf*ones(10,1);
    return
end

v = 1/2*x'*P_xi*x + 1/2*(x-xi)'*R*(x-xi);

if v < 0 
    y_dot = inf*ones(10,1);
    return
end


%% forwarding dynamics
y_dot = F*[x1;x2;x3;x4] + B*u;
xi_dot = -k * dv_dxi';
y_dot(5:8) = xi_dot;
y_dot(9) =  1/2 * x' * Q * x + 1/2 * u^2; 
y_dot(10) = - c_x_xi;


end

%% event to stop ode45
function [value, isterminal, direction] = myEvent(T, Y)    

% stop the code if nan appears, which means that the controls have not
% satisfied the inequality  

    if any(isnan(Y)) 
        condition = 0;
    else
        condition = 1;
    end

    value      = condition;
    isterminal = 1;   % Stop the integration
    direction  = 0;
    
end
