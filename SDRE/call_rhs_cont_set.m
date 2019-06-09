% Simulate the result of the SDRE control which uses the LQR methodolody at
% each time step

clear 
clc
close all

% parameters 
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

%% solver initial
tspan = [0,100];

% initial condtions
N0 = 1;
T0 = 0.2;
I0 = 0.15;
M0 = 0;     % this should zero at t = 0
U0 = 1;

% origin shift
x1_0 = N0 - 1/b2;
x2_0 = T0;
x3_0 = I0 - s/d1;
x4_0 = M0;
x5_0 = 0;

% input initial contditions into solver
y0 = [x1_0;x2_0;x3_0;x4_0;x5_0];


opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,y] = ode45(@(t,y) rhs_cont_set(t,y,1), tspan, y0,opts);

% obtain the control input
u = zeros(1,length(y));
for i = 1:length(y)
    [~,u(1,i)] = rhs_cont_set(t(i),y(i,:),1);
end

% shift back the origin
y(:,1) = y(:,1) + 1;
y(:,3) = y(:,3) + 1.65;


%% plots
fprintf('NO CONSTRAINT CASE \n')
fprintf('time taken to kill tumor = %f with cost = %f \n\n', min(t((y(:,2)<1e-4))),y(end,5))

plot(t,y(:,1:4),'-','linewidth',3)
hold on
plot(t,u,'-','linewidth',3)
plot(tspan,[0.75,0.75],'--k','linewidth',2)
legend('Normal Cell','Tumor Cell','Immune Cell','Drug Dose','Control Input','Normal Cell Limit')
grid on
xlabel('Time (days)')
ylabel('Normalised Cell Population')