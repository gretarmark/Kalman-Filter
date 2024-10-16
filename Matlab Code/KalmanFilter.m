%% Kalman Filter

clear all; close all; clc;

% Parameters
a1 = 1.2272;     %[cm2] Area of outlet pipe 1
a2 = 1.2272;     %[cm2] Area of outlet pipe 2
a3 = 1.2272;     %[cm2] Area of outlet pipe 3
a4 = 1.2272;     %[cm2] Area of outlet pipe 4

A1 = 380.1327;   %[cm2] Cross sectional area of tank 1
A2 = 380.1327;   %[cm2] Cross sectional area of tank 2
A3 = 380.1327;   %[cm2] Cross sectional area of tank 3
A4 = 380.1327;   %[cm2] Cross sectional area of tank 4

gamma1 = 0.45;  % Flow distribution constant. Valve 1
gamma2 = 0.40;  % Flow distribution constant. Valve 2

g = 981;        %[cm/s2] The acceleration of gravity
rho = 1.00;     %[g/cm3] Density of water

p = [a1; a2; a3; a4; A1; A2; A3; A4; gamma1; gamma2; g; rho];

% Simulation scenario
t0 = 0.0; % [s] Initial time
t_f = 20*60; % [s] Final time
Ts = 4; % [s] Sample time
t = [t0:Ts:t_f]'; % [s] Sample instants
N = length(t); 

m10 = 0.0; % [g] Liquid mass in tank 1 at time t0
m20 = 0.0; % [g] Liquid mass in tank 2 at time t0
m30 = 0.0; % [g] Liquid mass in tank 3 at time t0
m40 = 0.0; % [g] Liquid mass in tank 4 at time t0

F1 = 300; % [cm3/s] Flow rate from pump 1
F2 = 300; % [cm3/s] Flow rate from pump 2
F3 = 0.0; % Unknown stochastic variables (normal distributed) 
F4 = 0.0; % Unknown stochastic variables (normal distributed)

x0 = [m10; m20; m30; m40];

u = [F1 ; F2];
d = [F3(:,1) ; F4(:,1)];
u_rep = [repmat(F1,1,N); repmat(F2,1,N)];
d_rep = [repmat(F3,1,N); repmat(F4,1,N)];

a = p(1:4,1)';
A = p(5:8,1)';

xs0 = 5000*ones(4,1);   % Initial guess of steady state

xs = fsolve(@ModifiedFourTankSystemWrap, xs0, [], u , d , p) % Steady State
ys = FourTankSystemSensor(xs,p)
zs = FourTankSystemOutput(xs,p)

%% 
syms t u_1 u_2 x_1(t) x_2(t) x_3(t) x_4(t) F3 F4 d_1 d_2

x10 = xs(1);
x20 = xs(2);
x30 = xs(3);
x40 = xs(4);

Dx_1 = diff(x_1);
Dx_2 = diff(x_2);
Dx_3 = diff(x_3);
Dx_4 = diff(x_4);

eq1 = rho*gamma1*u_1 + rho*a3*sqrt(2*g*(x_3/(rho*A3))) - rho*a1*sqrt(2*g*(x_1/(rho*A1))) == Dx_1;
eq2 = rho*gamma2*u_2 + rho*a4*sqrt(2*g*(x_4/(rho*A4))) - rho*a2*sqrt(2*g*(x_2/(rho*A2))) == Dx_2;
eq3 = rho*(1 - gamma2)*u_2 - rho*a3*sqrt(2*g*(x_3/(rho*A3))) + rho*d_1 == Dx_3;
eq4 = rho*(1 - gamma1)*u_1 - rho*a4*(sqrt(2*g*(x_4/(rho*A4)))) + rho*d_2 == Dx_4;
[V,Y] = odeToVectorField([eq1;eq2;eq3;eq4;]);

V = strrep(string(V),'Y[1]',string(Y(1)));
V = strrep(string(V),'Y[2]',string(Y(2)));
V = strrep(string(V),'Y[3]',string(Y(3)));
V = strrep(string(V),'Y[4]',string(Y(4)));
    
V = [V(2) ; V(1) ; V(3) ; V(4)];
V = str2sym(V);
        
A = (jacobian(V,Y));
A = vpa([A(:,2) A(:,1) A(:,3) A(:,4)],8);
A = subs(subs(subs(subs(A,x_1,xs(1)),x_2,xs(2)),x_3,xs(3)),x_4,xs(4));
% A = subs(subs(subs(subs(A,x_1,ys(1)),x_2,ys(2)),x_3,ys(3)),x_4,ys(4));
A = vpa(A,7);
Ac = cast(A,'double')
    
B = (jacobian(V,[u_1 ; u_2]));
B = vpa(B,8);
Bc = cast(B,'double')
    
At = [A1 ; A2 ; A3 ; A4];
C=diag(1./(rho*At))
Cz=C(1:2,:)

E = jacobian(V,[d_1 ; d_2]);
E = vpa(E,4);
Ec = cast(E,'double')

%% 
% Qbar discretation, in the same place as c2dzoh e-a 
[Ad,Bd]=c2dzoh(Ac,Bc,Ts)
[~,Ed]=c2dzoh(Ac,Ec,Ts)

% noise covariance 
Qe = diag([5^2 ; 5^2]); 
Rv = diag([2^2;  2^2; 2^2; 2^2]);
% generate noise vector 
qe = chol(Qe,'lower')
W = qe*randn(2,N) % process noise
rv = chol(Rv,'lower')
V = rv*randn(4,N) % measurement noise 

% new noise W_bar 
Qe_bar = diag(Ed.^2*Qe*[1; 1])

% W_bar = zeros(4, N);
% for i =1:1:N
%     W_bar(:, i) = Ed*W(:, i); % process noise
% end

qe_bar = chol(Qe_bar,'lower')
W_bar = qe_bar*randn(4,N) % process noise

Gd_bar = diag([1;1;1;1]) % noise matrix  

T_plot = t0:Ts:t_f;


% stochastic simulation 
% new system :
% X_k+1 = Ad * x_k + Bd * u_k +  Ed*d_k + Gd_bar * W_bar_k;

% pre-allocate memory for storing results 
x_stochastic = zeros(4,N);
x_stochastic(:,1) = x0 - xs;
xbar_stochastic = x_stochastic; % simulation with W_bar

uss = [300 ; 300];
u_stochastic = u_rep - uss;
y_stochastic = zeros(4,N) - ys; 
ybar_stochastic = y_stochastic; % simulation with W_bar

% simulation for normal case 
for i=1:N-1
   x_stochastic(:,i+1) = Ad*x_stochastic(:,i) + Bd*u_stochastic(:,i) + Ed*(W(:,i)+d_rep(:,i));
   y_stochastic(:,i+1) = C*x_stochastic(:,i+1) + V(:,i+1);
   P = i;
end

for i=1:N-1
   xbar_stochastic(:,i+1) = Ad*xbar_stochastic(:,i) + Bd*u_stochastic(:,i) + Ed*d_rep(:,i) + Gd_bar*W_bar(:,i);
   ybar_stochastic(:,i+1) = C*xbar_stochastic(:,i+1) + V(:,i+1);
   P = i;
end

% plot 
close all;
figure;
sgtitle('Stochastic simulation (without KF )') 
subplot(2,1,1)
title("normal case")
plot(T_plot, y_stochastic' + ys', 'DisplayName', "normal stochatic simulation")
xlim([0 1200])

% We can see that the state values go to xs or the steady states.
% xs

subplot(2,1,2)
title("modified case")
plot(T_plot, ybar_stochastic' + ys', 'DisplayName', "modified stochatic simulation")
xlim([0 1200])

% you can notice that the modified stochastic case has almost the same
% performance as the normal case.

% static KF (if you want to test you can replace dynamic kf gain with static kalman filter gain )
% P_cov = dare(Ad', C', Qe_bar, Rv)
% R_e = C*P_cov*C' + Rv;
% fprintf('Kalman gain is:')
% K_f = P_cov*(C'/R_e)

% dynamic kf
P_d = diag([0.015^2, 0.015^2, 0.015^2, 0.015^2]) % initial guess on P 

% pre-allocate memory for storing results 
x_stochastic = zeros(4,N);
x_stochastic(:,1) = x0 - xs;
xbar_stochastic = x_stochastic; % simulation with W_bar

uss = [300 ; 300];
u_stochastic = u_rep - uss;
y_stochastic = zeros(4,N) - ys; % y measurement 
y_est = zeros(4, N)- ys;

ybar_stochastic = y_stochastic; % simulation with W_bar

% simulation for normal case 
for i=1:N-1
   P_d = Ad*P_d*Ad' + Gd_bar*Qe_bar*Gd_bar';
   K_df = P_d*C'*(C*P_d*C'+Rv); % kalman gain 
   % simulation 
   x_stochastic(:,i+1) = Ad*x_stochastic(:,i) + Bd*u_stochastic(:,i) + Ed*(W(:,i)+d_rep(:,i));
   y_stochastic(:,i+1) = C*x_stochastic(:,i+1) + V(:,i+1);
   e = y_stochastic(:,i+1) - C*x_stochastic(:,i+1);
   x_est = x_stochastic(:, i+1) + K_df*e;
   y_est(:, i+1) = C*x_est;
   P_d = P_d - K_df*(C*P_d*C'+Rv) *K_df';
end

% plot 
close all;
figure;
plot(T_plot, y_stochastic(1:2,1:end)' + ys(1:2)','.')
title('Stochastic simulation with dynamic Kalman filter') 
hold on 
plot(T_plot, y_est(1:2,1:end)' + ys(1:2)', 'LineWidth',2.5)
xlim([0 1200])
hold off 



%% Functions

function xdot = ModifiedFourTankSystem(t,x,u,d,p)

% FOURTANKSYSTEM Model dx/dt = f(t,x,u,p) for 4-tank System
%
% This function implements a differential equation model for the
% 4-tank system.
%
% Syntax: xdot = FourTankSystem(t,x,u,p)
% Unpack states, MVs, and parameters

m = x; % Mass of liquid in each tank [g]
F = [u; d]; % Flow rates in pumps [cm3/s]
a = p(1:4,1); % Pipe cross sectional areas [cm2]
A = p(5:8,1); % Tank cross sectional areas [cm2]
gamma = p(9:10,1); % Valve positions [-]
g = p(11,1); % Acceleration of gravity [cm/s2]
rho = p(12,1); % Density of water [g/cm3]

% Inflows
qin = zeros(4,1);
qin(1,1) = gamma(1)*F(1); % Inflow from valve 1 to tank 1 [cm3/s]
qin(2,1) = gamma(2)*F(2); % Inflow from valve 2 to tank 2 [cm3/s]
qin(3,1) = (1-gamma(2))*F(2); % + F(3); % Inflow from valve 2 to tank 3 [cm3/s]
qin(4,1) = (1-gamma(1))*F(1); % + F(4); % Inflow from valve 1 to tank 4 [cm3/s]

% Outflows
h = m./(rho*A); % Liquid level in each tank [cm]
qout = a.*sqrt(2*g*h); % Outflow from each tank [cm3/s]

% Differential equations
xdot = zeros(4,1);
xdot(1,1) = rho*(qin(1,1)+qout(3,1)-qout(1,1)); % Mass balance Tank 1
xdot(2,1) = rho*(qin(2,1)+qout(4,1)-qout(2,1)); % Mass balance Tank 2
xdot(3,1) = rho*(qin(3,1)-qout(3,1) + F(3)); % Mass balance Tank 3
xdot(4,1) = rho*(qin(4,1)-qout(4,1) + F(4)); % Mass balance Tank 4

end

function xdot = ModifiedFourTankSystemWrap(x,u,d,p)

xdot = ModifiedFourTankSystem(0,x,u,d,p);
end

function y = FourTankSystemSensor(x,p)

rho = p(12);
A = p(5:8);
y = zeros(4,1);
for i=1:4   
    y(i,1) = x(i,1)./(rho*A(i,1));
end
end

function y = FourTankSystemSensorUnit(x,p)

rho = p(12);
A = p(5:8);
y = zeros(4,1);
y(:,1) = x(:,1)./(rho*A(:,1));

end

function z = FourTankSystemOutput(x,p)

rho = p(12);
A = p(5:8);
z = zeros(2,1);
for i=1:2
    z(i,1) = x(i,1)./(rho*A(i,1));
end
end

function z = FourTankSystemOutputUnit(x,p)

rho = p(12);
A = p(5:8);
z = zeros(2,1);
z1 = x(1,1)/(rho*A(1,1));
z2 = x(2,1)/(rho*A(2,1));
z = [z1; z2]; 
end

function Qoutflow = FourTankSystemOutFlow(h,p)
% y = h = height
a = p(1:4);
g = p(11);
Qoutflow = zeros(4,1);
Qoutflow(:,1) = a.*sqrt(2*g*h(:,1));

end

function [Abar,Bbar]=c2dzoh(A,B,Ts)

[nx,nu]=size(B);
M = [A B; zeros(nu,nx) zeros(nu,nu)]*Ts;
Phi = expm(M);
Abar = Phi(1:nx,1:nx);
Bbar = Phi(1:nx,nx+1:nx+nu);

end

function [Abar,Qbar] = c2dsw(A,G,Ts)

nx = size(A,1);
nx1 = nx+1;
nx2 = nx+nx;
M = [-A G*G' ; zeros(nx,nx) A']*Ts;
phi = expm(M);
Abar = phi(nx1:nx2,nx1:nx2)';
Qbar = Abar*phi(1:nx,nx1:nx2);
end
