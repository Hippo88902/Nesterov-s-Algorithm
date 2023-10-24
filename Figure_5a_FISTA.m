clear 
clc
rng(100)
m = 100; % The number of data (row)
n = 500; % The number of variables (column)
A = normrnd(0, 1, [m,n]); % generate data
b = 0.01 * normrnd(0, 25, [m,1]); % generate random coefficient
la = 4;
K = 1000;
x0 = 0.1 * ones(n,1);
xk = x0;
xk_1 = x0;
yk = x0;
cost = zeros(K-1, 1);
L = max(eig(A' * A));
z = x0;
t1 = 1;
for i = 1: 10000
    z = (yk - 0.5 * L^-1 * A' * (A * yk - b));
    xk = max(abs(z) - la * L^-1, 0) .* sign(z);
    t2 = (1 + sqrt(1 + 4 * t1^2)) / 2;
    yk1 = xk + (t1 - 1) * t2^-1 * (xk - xk_1);
    cost(i) = 0.5 * norm(A * xk - b, 2)^2 + la * norm(xk, 1);
    t1 = t2;
    xk_1 = xk;
    yk = yk1;
end
f_opt = 0.5 * norm(A * xk - b, 2)^2 + la * norm(xk, 1);
uk = x0;
uk_1 = x0;
vk = x0;
z_ = x0;
t_1 = 1;
error = zeros(K-1,1);
for i = 1: K
    z_ = (vk - 0.5 * L^-1 * A' * (A * vk - b));
    uk = max(abs(z_) - la * L^-1, 0) .* sign(z_);
    t_2 = (1 + sqrt(1 + 4 * t_1^2)) / 2;
    vk1 = uk + (t_1 - 1) * t_2^-1 * (uk - uk_1);
    error(i) = abs(0.5 * norm(A * uk - b, 2)^2 + la * norm(uk, 1) - f_opt);
    t_1 = t_2;
    uk_1 = uk;
    vk = vk1;
end

iteration = linspace(1,K,K);
figure(1)
semilogy(iteration, error, 'color', 'red')
xlabel('iterations')
ylabel('error')
% A' * (b - A * x_k1)