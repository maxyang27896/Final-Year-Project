% This code is presents the result of a linear control using the LQR
% methodology, linearising the system about the origin. 

clear 
clc
close all

% tumor free equilibrium (desired)
[A,B,C,D] = cc_sys_linear(0,0,0,0,0);
lambda_1 = eig(A);

% solve the riccati equation for system about origin
Q = diag([0,1000,0,0.2]);
[X,L,G] = care(A,B,Q);

% calculate the K-gain required
K = -B'*X;


%% simulate the controlled system

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

tspan = [0,100];

% initial condtions, equilibrium point tumor free = [1,0,1.65];
N0 = 1;
T0 = 0.2;
I0 = 0.15;
M0 = 0;   

% origin shift
x1_0 = N0 - 1/b2;
x2_0 = T0;
x3_0 = I0 - s/d1;
x4_0 = M0;
x5_0 = 0;

% input initial parameters
y0 = [x1_0;x2_0;x3_0;x4_0;x5_0];


%% solving the differentiation with an input, run different tests
[t,y] = ode45(@(t,y) rhs_linear(t,y,K,Q), tspan, y0);

% shift back the origin
y(:,1) = y(:,1) + 1/b2;
y(:,3) = y(:,3) + s/d1;

plot(t,y(:,1:4),'-o')
grid on
ylabel('Number of Cells')
xlabel('Time (Number of days)')
legend('Normal Cell','Tumor Cell','Immune Cell','Drug Dose')

