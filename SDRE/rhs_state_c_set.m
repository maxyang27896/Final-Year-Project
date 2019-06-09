function [y_dot,u] = rhs_state_c_set(t,y,epsilon)
% define ODE dynamics here to be used in solver, states x = [N,I,T,M]' and
% control u = u(t)

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



% input y-vector 
x1 = y(1);
x2 = y(2);
x3 = y(3); 
x4 = y(4);

%% ==========================================================================================
% input SDRE LQR controller here, as a function of the states, one type of
% SDC parametrisation
A11 = -r2*(1+x1*b2);
A12 = -c4*(x1+1/b2);
A14 = -a3*(x1+1/b2);
A21 = -c3*x2;
A22 = r1*(1-b1*x2) - (c2*s/d1+c3/b2);
A23 = -c2*x2;
A24 = -a2*x2;
A32 = rho*(x3 + s/d1)/(alpha + x2) - c1*(x3 + s/d1) - x4;
A33 = -d1;
A34 = -a1*(x3 + s/d1) + x2;
A44 = -d2;

A = [A11,A12,0,A14;A21,A22,A23,A24;0,A32,A33,A34;0,0,0,A44];
B = [0;0;0;1];
C = eye(4,4);
D = zeros(4,1);


% input control using LQR 
sys = ss(A,B,C,D);


% add a state constraint to the penalties x1 > 0.75
N = 1;
phi_x1 = 1/(abs(x1 + 1/b2 - 0.75) + epsilon)^(2*N);
% phi_x4 = 1/(abs(x4 - 1)+ 0.001)^(2*N);
phi_x4 = 0;
delta_h = diag([-1,0,0,1]);


C = delta_h * A;
D = delta_h * B;
W = diag([phi_x1,0,0,phi_x4]);


% this is the additional penalty Q_c = C(x)'*W(x)*C(x)
Q_c =  C' * W * C;
Q = diag([0,100,0,0.1]) + Q_c;

% penalty on control 
R_c = D' * W * D;
R = 1.2 + R_c;


[K,S,e] = lqr(sys,Q,R);
u = -K * [x1,x2,x3,x4]';    

% % hard control bound 
% u = min(max(u,0),1);
 

%% forwarding the dyanmics 
y_dot = A * [x1;x2;x3;x4] + B * u;
y_dot(5) = 1/2 * [x1;x2;x3;x4]' * Q * [x1;x2;x3;x4] + 1/2*R*u^2;

end