function [y_dot] = rhs_linear(t,y,K,Q)

% input y-vector 
x1 = y(1);
x2 = y(2);
x3 = y(3); 
x4 = y(4);

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


%% ==========================================================================================
% input dynamic control law
u = -K*[x1;x2;x3;x4]; 


% non-linear dynamics
x1_dot = -r2*x1*(1+b2*x1) - c4*x2*x1 - c4/b2*x2 -a3*x4*x1 - a3/b2*x4;
x2_dot = r1*x2*(1-b1*x2) - (c2*s/d1 + c3/b2)*x2 - c3*x2*x1 - c2*x2*x3 - a2*x4*x2;
x3_dot =  - c1*s/d1*x2 - d1*x3 - a1*s/d1*x4 + rho*s/d1*x2/(alpha + x2) + rho*x2*x3/(alpha + x2) ...
    - c1*x3*x2 - a1*x4*x3;
x4_dot = u - d2*x4;


% defining the rates: y_dot = [x1_dot;x2_dot;x3_dot;x4_dot];
y_dot = [x1_dot;x2_dot;x3_dot;x4_dot];
y_dot(5) =  1/2*[x1;x2;x3;x4]'*Q*[x1;x2;x3;x4] + 1/2*u^2; % calculating the cost
end