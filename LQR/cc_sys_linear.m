function [ A,B,C,D ] = cc_sys_linear(x1_0,x2_0,x3_0,x4_0,U0)
% function to generate the linearised state space systems about points (N0,T0,I0,u0,v0)

syms x1 x2 x3 x4 u
% syms a1 a2 a3 b1 b2 alpha c1 c2 c3 c4 d1 d2 r1 r2 s rho
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
x1_dot = -r2*x1*(1+b2*x1) - c4*x2*x1 - c4/b2*x2 -a3*x4*x1 - a3/b2*x4;
x2_dot = r1*x2*(1-b1*x2) - (c2*s/d1 + c3/b2)*x2 - c3*x2*x1 - c2*x2*x3 - a2*x4*x2;
x3_dot =  - c1*s/d1*x2 - d1*x3 - a1*s/d1*x4 + rho*s/d1*x2/(alpha + x2) + rho*x2*x3/(alpha + x2) ...
    - c1*x3*x2 - a1*x4*x3;
x4_dot = u - d2*x4;

% jacobian 
A_syms(x1,x2,x3,x4,u) = [diff(x1_dot,x1), diff(x1_dot,x2), diff(x1_dot,x3),diff(x1_dot,x4);...
    diff(x2_dot,x1), diff(x2_dot,x2), diff(x2_dot,x3), diff(x2_dot,x4);...
    diff(x3_dot,x1), diff(x3_dot,x2), diff(x3_dot,x3), diff(x3_dot,x4);...
    diff(x4_dot,x1), diff(x4_dot,x2), diff(x4_dot,x3), diff(x4_dot,x4)];
B_syms(x1,x2,x3,x4,u) = [diff(x1_dot,u);diff(x2_dot,u);diff(x3_dot,u);diff(x4_dot,u)];


% state-space system
A=double(A_syms(x1_0,x2_0,x3_0,x4_0,U0));
B=double(B_syms(x1_0,x2_0,x3_0,x4_0,U0));
C=double(eye(4,4));
D=double(zeros(4,1));

end
