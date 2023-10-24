clear 
clc
rng(1)
m = 100; % The number of data (row)
n = 500; % The number of variables (column)
A = normrnd(0, 1, [m,n]); % generate data
b = normrnd(0, 25, [m,1]); % generate random coefficient
la = 4;
K = 20000;
x0 = 0 * ones(n,1);
x_k1 = x0;
y_k1 = x0;
cost = zeros(K-1,1);
L = max(eig(A' * A));
z = x0;
t1 = 1;

for i = 1: K
    z = (y_k1 - L^-1 * A' * (A * y_k1 - b));
    x_k2 = max(abs(z) - la * L^-1, 0) .* sign(z);
    t2 = (1 + sqrt(1 + 4 * t1^2)) / 2;
    y_k2 = x_k1 + (t1 - 1) * t2^-1 * (x_k2 - x_k1);
    cost(i) = 0.5 * norm(A * x_k1 - b, 2)^2 + la * norm(x_k1, 1);
    t1 = t2;
    x_k1 = x_k2;
    y_k1 = y_k2;
end

iteration = linspace(1,K,K);
figure(1)
semilogy(iteration, cost, 'color', 'red')
xlabel('iterations')
ylabel('cost')


% A' * (b - A * x_k1)
    
    