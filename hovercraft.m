% Builds and solves a simple linear program
clear all
close all

% rand('state',0);
% randn('state',0);

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
% cvx_optval: 3.7188

%% Test 3 - Optimal Control 
% minimize sum(f_cost(z)) + sum(f_cost(u))
% subject to: z(t+1) = A(t)*z(t) + B(t)*u(t), t = 0,..., N-1
% N: time horizon of the problem
% z: system state at time t (in dimension k)
% u: control action at time t (in dimension l)
% In this case, the objective is c'*x
% Optimization variable x for the LP: (u(0), z(1), u(1),...,u(N-1),z(N))

% Hovercraft:
A_state = [0 1 0 0;
           0 0 0 0;
           0 0 0 1;
           0 0 0 0];
       
B_state = [0 0;
           1 0;
           0 0;
           0 1];
C_state = [1, 1, 1, 1];
D_state = [0];
sys = ss(A_state, B_state, C_state, D_state);

z0 = [2;50;2;15];

k = 4; % Dimension of system state
l = 2; % Dimension of control action
N = 10; % Time Horizon
n = N*(k+l); % Length of optimization variable

Q = 10*eye(k);
R = 0.01*eye(l);
Qf = 10*eye(k);
euler = 1;

Qtilde = euler*Q;
Rtilde = euler*R;
Qftilde = Qf;

% Euler discretization:
% Atilde = eye(4) + euler*A_state;
% Btilde = euler*B_state;

% ZOH discretization: 
T = 0.01; % Sample time is 0.01 s
zoh_sys = c2d(sys, T, 'zoh');
Atilde = zoh_sys.A;
Btilde = zoh_sys.B;

%% Lateral Motion of a Car:
% V = 10; T = 0.01;
% Atilde = [1, 0; V*T, 1];
% Btilde = [T; 0.5*V*T^2];
% z0 = [-5; 0];
% 
% k = 2; % Dimension of system state
% l = 1; % Dimension of control action
% N = 100; % Time Horizon
% n = N*(k+l); % Length of optimization variable
% 
% Q = 1*eye(k);
% R = 1*eye(l);
% Qf = 1*eye(k);
% euler = 1;
% 
% Qtilde = euler*Q;
% Rtilde = euler*R;
% Qftilde = Qf;

% Quadratic Objective

cvx_begin quiet

    variables x(k, (N+1)) u(l, N);
    expressions s_x(N) s_u(N) s_f;
    for i = 1:N
        s_x(i) = x(:,i)'*Qtilde*x(:,i);
        s_u(i) = u(:,i)'*Rtilde*u(:,i);
    end
    s_f = x(:,N+1)'*Qftilde*x(:,N+1);
    minimize sum(s_x) + sum(s_u) + s_f + (z0'*Qtilde*z0)
    subject to
    for i = 1:N
        x(:,i+1) == Atilde*x(:,i) + Btilde*u(:,i);
    end

    % Initial conditions
    u(1,:) <= 2500;
    u(2,:) <= 2500;
    -u(1,:) <= 2500;
    -u(2,:) <= 2500;
    x(:,1) == z0;
    x(:,N+1) == [0;0;0;0];
    

cvx_end

% Linear Objective
% cvx_begin quiet
%    variables x(k, (N+1)) u(l, N);
%    expressions s_x(N) s_u(N) s_f;
%    for i = 1:N
%        s_x(i) = norm(Qtilde*x(:,i),1);
%        s_u(i) = norm(Rtilde*u(:,i),1);
%    end
%    s_f = norm(Qf*x(:,N+1),1);
%    minimize sum(s_x) + sum(s_u) + s_f
%    subject to
%       for i = 1:N
%            x(:,i+1) == Atilde*x(:,i) + Btilde*u(:,i);
%       end
%    x(:,1) == z0; % x(0) = z(0)
%    x(:,N+1) == [0;0;0;0];
%    
% cvx_end

% Extract control inputs and states:
u = zeros(l, N);
x = zeros(k, N+1);
total_penalty = cvx_optval;
penalty = [sum(s_x); sum(s_u); sum(s_f)];


for i = 1:N
    x(:,i) = cvx_optpnt.x(:,i);
    u(:,i) = cvx_optpnt.u(:,i); 
end
x(:, N+1) = cvx_optpnt.x(:,N+1);

figure(1)
title('Hovercraft trajectory');
xlabel('X');
ylabel('Y');
grid on
hold on
dest = plot(0,0, 'r*');
start = plot(z0(1), z0(3), 'b^');
plot(z0(1), z0(3), 'k-');
set(dest,'linewidth',2);
set(start,'linewidth',1.5);
for i = 1:N+1
    pause(0.25);
    start = plot(x(1,i), x(3,i), 'b^');
    set(start,'linewidth',1.5);
end

% Trajectory in direction of motion of car:
% x_init = 2;
% xtr = zeros(size(x(2,i)));
% for i = 1:N+1
%     xtr(i) = x_init + V*T*i;
% end
% 
% figure(1)
% title('Car trajectory');
% xlabel('X');
% ylabel('Y');
% grid on
% hold on
% % dest = plot(0,0, 'r*');
% start = plot(x_init, z0(2), 'b^');
% % set(dest,'linewidth',2);
% set(start,'linewidth',1.5);
% for i = 1:N+1
%     pause(0.25);
%     start = plot(xtr(i), x(2,i), 'b^');
%     set(start,'linewidth',1.5);
% end


%% Linear Objective

% cvx_begin quiet
%    variables x(k, (N+1)) u(l, N);
%    expressions s_x(N) s_u(N) s_f;
%    for i = 1:N
%        % s_x(i) = x(: ,1+(i-1) : i*k)'*Q*x(1+(i-1)*k : i*k);
%        % s_u(i) = u(1+(i-1)*l : i*l)'*R*u(1+(i-1)*l : i*l);
%        s_x(i) = norm(Qtilde*x(:,i),1);
%        s_u(i) = norm(Rtilde*u(:,i),1);
%    end
%    s_f = norm(Qf*x(:,N+1),1);
%    minimize sum(s_x) + sum(s_u) + s_f
%    subject to
%       for i = 1:N
%            x(:,i+1) == Atilde*x(:,i) + Btilde*u(:,i);
%       end
%    x(:,1) == z0; % x(0) = z(0)
%    x(:,N+1) == [0;0;10;0];
%    u(:,3) == [1;0];
%       % x(1+N*k : (N+1)*k) == [0; 0; 0; 0]; % Final state is at origin.
% cvx_end


% Setting up feasibility constraint for an LP
% row = N*k;
% col = N*(k+l) + k; % See pg. 552 BV
% 
% A = zeros(row,col);
% b = zeros(row,1);
% set = [-A_state, -B_state, eye(k,k)];
% len = length(set);
% for i = 1:N
%     A(1+(i-1)*k:i*k ,1+(i-1)*(k+l):len + (i-1)*(k+l)) = set;
% end
% 
% A = A(:,k+1:end);
% b(1:k) = A_state*z_0;
% c = ones(n,1);