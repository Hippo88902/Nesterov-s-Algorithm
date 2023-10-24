clc;clear;close all;
rng(1)
m = 100; % The number of data (row)
n = 500; % The number of variables (column)
A = normrnd(0, 1, [m,n]); % generate data
b = normrnd(0, 25, [m,1]); % generate random coefficient
L = max(eig(A' * A));
la = 4;
K = 15000;
x0 = 0.5 * ones(n,1);
xk = x0;
xk_1 = x0;
yk = x0;
yk_1 = x0;
cost = zeros(K-1,1);
error = zeros(K-1,1);
r = 4;
s = 1 / L ;

for j = 1: 10000
    beta = (j-1)/(j+r-1);
    u = yk_1 - s * A' *(A * yk_1 - b);
    argmin = max( abs(u) - la*s, 0 ) .* sign(u);
    Gs = (yk_1 - argmin) / s;
    xk = yk_1 - s * Gs ;
    yk = xk + beta * (xk - xk_1);
    cost(j) = 0.5 * norm(A * xk - b)^2 + la * norm(xk, 1);
    xk_1 = xk;  
    yk_1 = yk;
end

uk = x0;
uk_1 = x0;
vk = x0;
vk_1 = x0;
f_opt = 0.5 * norm(A * xk - b)^2 + la * norm(xk, 1);
for i = 1: K
    beta = (i-1)/(i+r-1);
    u = vk_1 - s * A' *(A * vk_1 - b);
    argmin = max( abs(u) - la * s, 0 ) .* sign(u);
    Gs = (vk_1 - argmin) / s;
    uk = vk_1 - s * Gs ;
    vk = uk + beta * (uk - uk_1);
    error(i) = abs(0.5 * norm(A * uk - b)^2 + la * norm(uk, 1) - f_opt);
    uk_1 = uk;  
    vk_1 = vk;
end

iteration = linspace(1,K,K);

figure(1)
semilogy(iteration,error,'color','red')
xlabel('iterations')
ylabel('f-f*')

