clear,close all,clc
% A comparision of NNLS using various iterative first order methods
%{  
% Descriptions :
% NNLS model 
%      min_{ x \geq 0} (1/2) || Ax - b ||_2^2
%
% Equivalent NNQP model
%      min_{ x \geq 0} (1/2) x'Qx - p'x 
% where Q = A'*A, p = A'*b
%
% For THETA and Theta,
%
% THETA = I_n - Q/norm(Q,2)
% Theta = p/norm(Q,2)
%
% Algorithms :
% PGD : projected gradient descent
% APGD : accelerated gradient descent, with variations :
%     1. Nesterov's parameter, and with restart
%     2. Paul Tseng's parameter, and adaptive restart
%     3. Constant Beta using conditional number, and adaptive restart 
% All the gradient descent use step size 1/norm(Q) (i.e.constant step size)
% It is possible to modify the code to extact line search step size
%
% Also see https://angms.science/doc/NMF/nnls_pgd.pdf
% Created : 2018-11-3, 
% Last update : 2018-11-7
% Andersen Ang @ UMONS, BE, angms.science
%}
%%
% dimensions
m = 100;  % row 
n = 10;   % column
% random generate matrix A
A = rand(m,n);
% random generate the ground true nonnegative vector x
xtrue = rand(n,1);
xtrue(xtrue<0) = 0;
% form the vector b
b = A * xtrue;
%% Some useful variables
% Initial guess
 x0 = rand(n,1);
% iteration max
 kmax = 100;
%% PGD 
 % Parameter matrix and vector
 THETA = eye(n) - A'*A/norm(A'*A,2);
 theta = A'*b/norm(A'*A,2);
 
 x = x0;  % initialize
 for k = 1 : kmax  % Loop
    x = max(THETA*x + theta, 0);  % update
    e_pgd(k) = norm(A*x - b);     % compute error
 end
%% APGD-1
 % Acceleration parameter
 alpha = 0.9; % alpha in (0, 1)
 % initialize the variable and the second variable 
 x = x0;
 y = x;
 % Loop
 for k = 1 : kmax
    % store previous variable
    x_ = x; 
    % store previous alpha
    alpha_ = alpha;
    % computation of acceleration parameter
    alpha = 0.5*(sqrt(alpha^4+4*alpha^2) - alpha^2);
    beta = alpha_*(1-alpha_)/(alpha_^2+alpha);
    % update
    x = max(THETA*y + theta, 0);
    % extrapolation
    y = x + beta*(x - x_);
    % compute error
    e_apgd(k) = norm(A*x - b);
 end
%% APGD-1 with restart
 % Acceleration parameter
 alpha = 0.9; % alpha in (0, 1)
  % store initial alpha value for restart
 alpha0 = alpha; 
 % initialize the variable and the second variable 
 x = x0;
 y = x;
 % Loop
for k = 1 : kmax 
    % store previous variable
    x_ = x;
    % store previous alpha
    alpha_ = alpha;
    % computation of acceleration parameter
    alpha = 0.5*(sqrt(alpha^4+4*alpha^2) - alpha^2);
    beta = alpha_*(1-alpha_)/(alpha_^2+alpha);
    % update
    x = max(THETA*y + theta, 0);
    % extrapolation
    y = x + beta*(x - x_);
    % compute error
    e_apgd_r(k) = norm(A*x - b);
    % check does error increase
    if (k > 1) && (e_apgd_r(k) > e_apgd_r(k-1))
      % restart, just consider a gradient step update
      x  = max(THETA*x_ + theta, 0);
      e_apgd_r(k) = norm(A*x - b);
      % restart 
      y = x;
      alpha = alpha0;
    end 
end
%% Plot errors vs iterations
semilogy(e_pgd,'k--','linewidth',2),hold on,
semilogy(e_apgd,'b','linewidth',6),hold on,
semilogy(e_apgd_r,'c','linewidth',2),hold off
legend('normal PGD','A','A + r');
title('NNLS','interpreter','latex','fontsize',20)
xlabel('Iteration/$k$','interpreter','latex','fontsize',15)
ylabel('$\| Ax - b\|$','interpreter','latex','fontsize',15)
axis tight,grid on