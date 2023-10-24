clear,close all,clc
%%
rng(1000)
% dimensions
m = 100;  % row 
n = 500;  % column
% random generate matrix A
A = normrnd(0, 1, [m,n]);
% random generate the ground true nonnegative vector x
xtrue = rand(n,1);
xtrue(xtrue<0) = 0;
% form the vector b
b = A * xtrue;
 %% Some useful variables
% Initial guess
 x0 = rand(n,1);
% iteration max
 kmax = 500;
% Parameter matrix and vector
 THETA = eye(n) - A' * A/norm(A'*A,2);
 theta = A' * b/norm(A'*A,2);
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
    alpha = 0.5*(sqrt(alpha^4 + 4*alpha^2) - alpha^2);
    beta = alpha_*(1 - alpha_)/(alpha_^2 + alpha);
    % update
    x = max(THETA*y + theta, 0);
    % extrapolation
    y = x + beta*(x - x_);
    % compute error
    e_apgd(k) = norm(A*x - b);
 end
 
 semilogy(e_apgd,'b'),hold on,
 legend('normal PGD','A','A + r');
 title('NNLS','interpreter','latex','fontsize',20)
 xlabel('Iteration/$k$','interpreter','latex','fontsize',15)
 ylabel('$\| Ax - b\|$','interpreter','latex','fontsize',15)
 axis tight,grid on

