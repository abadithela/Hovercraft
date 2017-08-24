function [xopt, lam2hist] = centerstep(A, b, c ,x0, t, alpha, beta, EPSTOL)
% -------- INPUTS -------
% x0 - Starting point for Newton's method
% A, b, c - Matrices for this problem
% t - Strength of indicator function in LP formulation
% alpha, beta - Parameters for back-tracking line search
% EPSTOL - Tolerance for exiting centering step method
% -------- OUTPUTS ------
% xopt - Solution from Newton's method
% count - Number of Newton iterations
% TODO: nuopt - not useful, need to remove

% TOLERANCES
FPTOL = 1e-3;
MAX = 40; 
xopt = []; lam2hist = [];

% Asserting feasibility:
% assert(min(x0) > 0);
% assert(norm(A*x0 - b) < FPTOL)
    

% Setup: 
count = 0;
x = x0;
c = t*c;
[row, col] = size(A);
% row = length(b)
% col = length(x0)

while(count  < MAX)
    
	%assert(min(x) > FPTOL);
    %assert(norm(A*x - b) < FPTOL);
    
    count = count  + 1;
    % Hessian and gradient:
    H = diag(1./x.^2, 0);
    g = c - 1./x;
    
    % Direct solve KKT to find search direction:
    kkt =  [H, A'; A, zeros(row, row)];
    NT = kkt\[-g; zeros(row,1)];
    xnt = NT(1:col);
    
    % Block Elimination
%     h = zeros(row, 1);
%     h1 = H\A';
%     h2 = H\g;
%     S = -A*h1; % Schur Complement
%     w = S\(A*h2 - h);
%     xnt = H\(-A'*w - g);
    
    % Newton decrement squared:
    lam2 = -g'*xnt;
    lam2hist = [lam2hist lam2/2];
    
    % assert(lam2 > FPTOL);
    % Exit from Newton's method
    if(lam2/2 < EPSTOL)
        break;
    end
    
    % Line Search
    lst = 1;
    
    while(min(x + lst*xnt) <= 0)
        lst = lst * beta;        
    end
        
    while(c'*lst*xnt - sum(log(x + lst*xnt)) + sum(log(x)) - alpha*lst*g'*xnt > FPTOL)
        lst = lst * beta;
    end
    
    % Update:
    x = x + lst*xnt;   
    
end

xopt = x;

end




























% fptol = 1e-10;
% assert(min(x0) > fptol);
% assert(norm(A*x0 - b) < fptol);
% count = 0; % Keeps track of number of Newton steps required for each iteration
% n = length(x0);
% f = @(x) t*c'*x - sum(log(x)); % Objective: tc'x + phi(x)
% 
% x = x0;
% c = t*c; % Multiply c by a factor of t for every outer iteration
% while(count<50)
%     count = count + 1;
%     H = diag(1./x.^2, 0);
%     g = c - 1./x;
%     D = diag(x.^2, 0);
%     nunt = (A*D*A')\(-t*A*D*c + b);
%     xnt = -t*D*c + x - D*A'*nunt;
%     lam2 = -xnt'*g;
%     
%     % Stopping Criterion for Newton iterations:
%     if(lam2/2 <= epsilon)
%         break;
%     end
%     
%     % Backtracking line search:
%     ls_t = 1; % ls_t = line search t. Setting line search parameter t to 1
%     while((x + ls_t*xnt) <= 0) % Making sure that x is in the domain of f. Since we have a log function, we require x+ls_t*xnt > 0
%         ls_t = ls_t*beta;
%     end
%     
%     while(f(x + ls_t*xnt) > f(x) + alpha*ls_t*g'*xnt)
%         ls_t = ls_t*beta;
%     end
%     
%     % Updating Newton Step:
%     x = x + ls_t*xnt;
% end
% 
% xopt = x; % Optimal x for current outer iteration for the given t.
% nuopt = nunt; % Optimal dual nu for the current outer iteration
% 
% end