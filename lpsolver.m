function [x_star, optval] = lpsolver(A, b, c)
%% ======== Standard Form LP Solver using Interior Point Methods ======== %
% Apurva Badithela (6/12/17) 
%
% minimize c'*x
% subject to: Ax = b, x>= 0
%
% Variables: x
% Implements Phase 1 and Phase 2 to find optimal x
% Input: A, b, c
% Output: Optimal point (x_star) and optimal value of objective (c'*x_star)
% Figures: Duality Gap vs. Number of Newton iterations

%% ===== PARAMETERS AND TOLERANCES ===== %
ALPHA = 0.03;
BETA = 0.5;
EPSTOL = 1e-6;
TTOL = 1e-3;
[m, n] = size(A);

%% ===== SETUP ===== %

% rand('state',0);
% randn('state',0);
% 
% % % Feasible:
% % A = [randn(m-1,n); ones(1,n)];
% % x0 = rand(n,1) + 0.1;
% % b = A*x0;
% % c = rand(n,1);
% 
% % Possibly Infeasible: 
% % A = [1 -2 -2  4; 
% %      5  3  4 -9;
% %      2  5 -3  6];
% % b = [-5; 14; 9];
% % c = rand(4,1);
% % 
% % [m,n] = size(A);
% f = @(x) c'*x;

x_star = [];
optval = [];
%% ===== Phase 1 ===== %

% Feasible x0 for Phase 1
x0 = A\b;

% Implementing Phase 1
% If min (x0) < 0, then need to implement Phase 1
% If min (x0) > 0, then skip Phase 1 and proceed to Phase 2

% Change of coordinates: z = x + (t-1)*ones(n,1)
if(min(x0) <= 0)
    t0 = 2 -min(x0);
    z0 = x0  + (t0-1)*ones(n,1);
    z_init = [z0;t0];
    
    % minimize c'*z_star
    % subject to: [A -A*ones(n,1)] z_star = b - A*ones(n,1)
    % z_star >= 0 , where z_star = [x; t]
    
    A_ph1 = [A, -A*ones(n,1)];
    b_ph1 = b - A*ones(n,1);
    c_ph1 = [zeros(n,1); 1];
    
    mu_ph1 = 20;
    t_ph1 = 1;
    gap_ph1 = (n+1)/t_ph1;
    gap1 = []; history_ph1 = [];
    
    while(gap_ph1 > EPSTOL)
        [z_ph1,lam2_ph1] = centerstep(A_ph1, b_ph1, c_ph1, z_init, t_ph1, ALPHA, BETA, EPSTOL);
        if(isempty(z_ph1))
            sprintf('Problem is infeasible.');
            return;
        end
        
        z_init = z_ph1;
        t_ph1 = mu_ph1*t_ph1;
        gap_ph1 = (n+1)/t_ph1;
        history_ph1 = [history_ph1 length(lam2_ph1)]; gap1 = [gap1 gap_ph1];
    end
    
    % Problem is: feasible if t < 1 at end of Phase 1
    %             infeasible, otherwise
    
    if(z_ph1(end) - 1 < TTOL)
        sprintf('Problem is feasible. Implementing Phase 2 ... ');
    else
        sprintf('Problem is infeasible.');
        return;
    end
    x_init = z_ph1(1:n) - (z_ph1(n+1)-1)*ones(n,1); 
% assert(z_ph1(end) - 1 < TTOL);

% Initial x for Phase 2
else
    x_init = x0;
end
%% ===== Phase 2 ===== %

t_ph2 = 1;
mu_ph2 = 20;
gap_ph2 = n/t_ph2;
hist_ph2 = [];
gap2 = [];

while(gap_ph2 > EPSTOL)
    [x_ph2, lam2_ph2] = centerstep(A, b, c, x_init,t_ph2, ALPHA, BETA, EPSTOL);
    if(isempty(x_ph2))
        sprintf('Problem is infeasible.');
        return;
    end
    x_init = x_ph2;
    t_ph2 = mu_ph2*t_ph2;
    gap_ph2 = n/t_ph2;
    hist_ph2 = [hist_ph2 length(lam2_ph2)]; gap2 = [gap2 gap_ph2];
end

figure (1)
[xx, yy] = stairs(cumsum(hist_ph2),gap2);
semilogy(xx,yy,'bo-');
xlabel('Newton Iterations');
ylabel('Duality Gap');
title('LP Solver - Newton iterartions per step');

%% ===== Verify with CVX ===== %

% Objective
f = @(x) c'*x;
solver_optval = f(x_init);
sprintf('Optimal value found by this solver: %f \n',f(x_init));

% cvx_begin
% variable x(n)
%    dual variables y z
%    minimize(c' * x)
%    subject to
%       y : A * x == b;
%       z : x >= 0;
% cvx_end
% 
% sprintf('Optimal value found by CVX: %f \n',cvx_optval);

x_star = x_init;
optval = solver_optval;

end