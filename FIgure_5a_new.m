clc;clear;close all;

rng(1)
m = 100;  % The number of data (row)
n = 500;  % The number of variables (column)
A = normrnd(0,1,[m,n]); % generate data
b = normrnd(0,25,[m,1]);

%x_1 = (A'*A)^(-1)*A'*b;
x_0 = zeros(n,1);
x_0(1) = 25;
%% FISTA
la = 4;
K = 10000;
L1 = max(eig(A' * A));
x0 = randn(n,1);
x_k1 = x0;
y_k1 = x0;
f_fista = zeros(K,1);
err_fista = zeros(K,1);
z = x0;
t1 = 1;
fopt = 1/2 * norm(A * x_0 - b)^2 + la * norm(x_0,1);

for i = 1: K
    z = (y_k1 - L1^-1 * A' * (A * y_k1 - b));
    x_k2 = max(abs(z) - la * L1^-1, 0) .* sign(z);
    t2 = (1 + sqrt(1 + 4 * t1^2)) / 2;
    y_k2 = x_k1 + (t1 - 1) * t2^-1 * (x_k2 - x_k1);
    f_fista(i) = 0.5 * norm(A * x_k1 - b, 2)^2 + la * norm(x_k1, 1);
    err_fista(i) = f_fista(i) - fopt;
    t1 = t2;
    x_k1 = x_k2;
    y_k1 = y_k2;
end

%% ISTA
rng(1)
L2 = norm(A)^2;
t = 1/L2;
%
f_ista = zeros(K,1);
err_ista = zeros(K,1);

for i = 1:K
    x0 = max(abs(x0 - t * A' * (A * x0-b))-la * t,0).*sign(x0-t*A'*(A*x0-b));
    f_ista(i) = 0.5 * norm(A * x0 - b) ^ 2 + la * norm(x0, 1);
    err_ista(i) = f_ista(i) - fopt;   
end
%% graph
semilogy(err_ista,'r')
%legend('FISTA','ISTA')
hold off
%x0
%A'*(A*x0-b)