%% ====== Examples using LP solver ====== %
% Apurva Badithela - 6/12/17
% Initialize A, b, and c matrices to be input to lpsolve.m
clear all
close all

rand('state',0);
randn('state',0);


%% Real World Problem 1 - Hovercraft Dynamics
% Hovercraft moving in the x,y plane with thrusts T_x and T_y controlling
% the motion of the hovercraft in the x and y directions:

% n = 32;
% A_state = [0 1 0 0;
%            0 0 0 0;
%            0 0 0 1;
%            0 0 0 0];
% B_state = [0 0;
%            1 0;
%            0 0;
%            0 1];
%        
% % Controllability Matrix:
% 
% Ctrl = [B_state, A_state*B_state, A_state*A_state*B_state, A_state*A_state*A_state*B_state];
% H = [A_state*A_state*A_state*B_state, A_state*A_state*B_state, A_state*B_state, B_state];
% x_des = [0; 0; 0; 0];
% 
% % A, b, c matrices:
% A = [H/2  -H/2  zeros(4,8) zeros(4,8)];
% b = x_des;
% c = [0.5*ones(8,1); 0.5*ones(8,1); zeros(8,1); zeros(8,1)];
%% Test 1
% clear all
% m = 100;
% n = 500;
% A = [randn(m-1,n); ones(1,n)];
% v = rand(n,1) + 0.1;
% b = A*v;
% c = randn(n,1);

%% Test 2
% clear all
% m = 3;
% n = 5;
% A = [-1 2 2 5 7;
%      8 3 -6 -9 4;
%      2 3 4 5 1];
% b = [22; -15; 24];
% c = rand(n,1);
% opt_sol = 3.7188

%% Linear MPC solved using LP:
mass = 1; % kg

% State: [x, y, x_dot, y_dot]
A_state = [zeros(2,2), eye(2); zeros(2,2), zeros(2,2)];        
B_state = [zeros(2,2); 1.0/mass*eye(2)];       
C_state = [1, 1, 1, 1];
D_state = [0];
sys = ss(A_state, B_state, C_state, D_state);
z0 = [2;2;50;15];
zf = [0; 0; 0;0];

k = 4; % Dimension of system state
l = 2; % Dimension of control action
N = 10; % Time Horizon

% Weights are constant for the entire time horizon
Q = 10*zeros(k);
R = 0.01*eye(l);
Qf = 5*eye(k);
euler = 1;

Qtilde = sqrtm(Q);
Rtilde = sqrtm(R);
Qftilde = sqrtm(Qf);

% Euler discretization:
% Atilde = eye(4) + euler*A_state;
% Btilde = euler*B_state;

% ZOH discretization: 
T = 0.01; % Sample time is 0.01 s
zoh_sys = c2d(sys, T, 'zoh');
Atilde = zoh_sys.A;
Btilde = zoh_sys.B;

row = N*k;
col = N*(k+l);
t = N*(k+l) + k; % See pg. 552 BV

A_dyn = zeros(row,t);
b_dyn = zeros(row,1);

set = [-Atilde, -Btilde, eye(k,k)];
len = length(set);
for i = 1:N
    A_dyn(1+(i-1)*k:i*k,1+(i-1)*(k+l):len+(i-1)*(k+l)) = set;
end

A_dyn = A_dyn(:,k+1:end);
b_dyn(1:k) = Atilde*z0;
c_dyn = zeros(k, col);
c_dyn(:, N*(k+l)-k+1:N*(k+l)) = eye(k);

% Finding Qftilde
Qftilde_dyn = zeros(k, col);
Qftilde_dyn(:, col-k+1:col) = Qftilde;
I_Qf = zeros(k, 2*t);
s1_p = zeros(k, 2*t);
s2_p = zeros(k, 2*t);
I_Qf(:, 1:k) = -eye(k);
I_Qf(:, k+1:2*k) = eye(k);
s1_p(:, 1:k) = eye(k);
s2_p(:, k+1:2*k) = eye(k);

% Finding Q's
Q_i = zeros(k, col, N);
I_Qi = zeros(k, 2*t, N);
s1_Qi = zeros(k, 2*t, N);
s2_Qi = zeros(k, 2*t, N);

for i = 1:N-1
    Q_i(:,(k+l)*(i-1)+l+1 : (k+l)*i,i) = Qtilde;
    I_Qi(:, 2*k+(i-1)*k +1 : 2*k + i*k ,i) = -eye(k);
    I_Qi(:, N*k + 2*k+(i-1)*k +1 : N*k + 2*k + i*k ,i) = eye(k);
    s1_Qi(:, 2*k+(i-1)*k+1 : 2*k + i*k ,i) =  eye(k);
    s2_Qi(:, N*k + 2*k+(i-1)*k+1 : N*k + 2*k + i*k ,i) =  eye(k);
end

% Finding R's
R_i = zeros(l, col, N);
I_Ri = zeros(l, 2*t, N);
s1_Ri = zeros(l, 2*t, N);
s2_Ri = zeros(l, 2*t, N);

for i = 1:N
    R_i(:,(k+l)*(i-1)+1 : (k+l)*(i-1) + l,i) = Rtilde;
    I_Ri(:,2*k+2*N*k +(i-1)*l+1:2*k+2*N*k+i*l, i) = -eye(l);
    I_Ri(:,2*k+2*N*k + N*l +(i-1)*l+1:2*k+N*l+2*N*k+i*l,i) = eye(l);
    s1_Ri(:,2*k+2*N*k+(i-1)*l+1 : 2*k+2*N*k+i*l,i) = eye(l);
    s2_Ri(:,2*k+2*N*k+N*l+(i-1)*l+1 : 2*k+2*N*k+N*l+i*l,i) = eye(l);
end

% Making bands in A:
A = [Qftilde_dyn, -Qftilde_dyn, I_Qf, s1_p;
     -Qftilde_dyn, Qftilde_dyn, I_Qf, s2_p];
for i = 1:N-1
    band = [Q_i(:,:,i), -Q_i(:,:,i), I_Qi(:,:,i), s1_Qi(:,:,i);
            -Q_i(:,:,i), Q_i(:,:,i), I_Qi(:,:,i), s2_Qi(:,:,i)];
    A = [A; band];
end

for i = 1:N
    band = [R_i(:,:,i), -R_i(:,:,i), I_Ri(:,:,i), s1_Ri(:,:,i);
            -R_i(:,:,i), R_i(:,:,i), I_Ri(:,:,i), s2_Ri(:,:,i)];
    A = [A; band];
end

A_band = [A_dyn, -A_dyn, zeros(row, 2*(2*k + 2*N*(k+l)))];
c_band = [c_dyn, -c_dyn, zeros(k, 2*(2*k + 2*N*(k+l)))];
A = [A; A_band; c_band];

% Setting up b:
[ROW, COL] = size(A);
b = zeros(ROW,1);
b(ROW - k - row + 1: ROW-k) = b_dyn;
b(ROW-k+1:ROW) = zf;

% Setting up c:
c = zeros(COL,1);
c(2*N*(k+l)+1:4*N*(k+l)+2*k) = ones(2*k+2*N*(k+l),1);


%% CVX Check

cvx_begin quiet
    variable x(COL,1);
    minimize (c'*x)
    subject to
    A*x == b
    -x <= 0
cvx_end

% Extracting inputs and states:
% x = x(1:col);
% u = zeros(l, N);
% z = zeros(k, N);
% 
% for i = 1:N
%     loc = 1 + (i-1)*(k+l);
%     u(:,i) = x(loc : loc + l-1);
%     z(:,i) = x(loc+l : loc+l + k - 1); 
% end
xu = x(1:col) - x(col+1:2*col);
u = zeros(l, N);
z = zeros(k, N+1);
for i = 1:N
    loc = 1 + (i-1)*(k+l);
    u(:,i) = xu(loc : loc + l-1);
    z(:,i) = xu(loc+l : loc+l + k - 1); 
end
z(:,N+1) = zf;

figure(100)
title('Hovercraft trajectory');
xlabel('X');
ylabel('Y');
grid on
hold on
dest = plot(0,0, 'r*');
start = plot(z0(1), z0(2), 'b^');
plot(z0(1), z0(2), 'k-');
for i = 1:N+1
    pause(0.25);
    start = plot(z(1,i), z(2,i), 'b^');
end

%% Test LP Solver

[x_sol, opt_sol] = lpsolver(A,b,c);
% Extracting data from LP solver:
xu_sol = x_sol(1:col) - x(col+1:2*col);
u_sol = zeros(l, N);
z_sol = zeros(k, N+1);
for i = 1:N
    loc = 1 + (i-1)*(k+l);
    u_sol(:,i) = xu_sol(loc : loc + l-1);
    z_sol(:,i) = xu_sol(loc+l : loc+l + k - 1); 
end
z_sol(:,N+1) = zf;

% figure(200)
% title('LP Solver: Hovercraft trajectory');
% xlabel('X');
% ylabel('Y');
% grid on
% hold on
% dest = plot(0,0, 'r*');
% start = plot(z0(1), z0(2), 'b^');
% plot(z0(1), z0(2), 'k-');
% for i = 1:N+1
%     pause(0.25);
%     start = plot(z_sol(1,i), z_sol(2,i), 'b^');
% end


% Extract inputs and states:
% u = zeros(l, N);
% z = zeros(k, N);
% optval = opt_sol;
% 
% for i = 1:N
%     loc = 1 + (i-1)*(k+l);
%     u(:,i) = x_sol(loc : loc + l-1);
%     z(:,i) = x_sol(loc+l : loc+l + k - 1); 
% end