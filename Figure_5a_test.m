clear 
clc
rng(1)
m = 100; % The number of data (row)
n = 500; % The number of variables (column)
A = normrnd(0, 1, [m,n]); % generate data
L = max(eig(A' * A));
A = A / sqrt(L);
b = normrnd(0, 25, [m,1]); % generate random coefficient
b = b / sqrt(L);
L = max(eig(A' * A));
la = 4;
K = 5000;
x0 = 1 * ones(n,1);
x_k1 = x0;
y_k1 = x0;
cost1 = zeros(K-1,1);
cost2 = zeros(K-1,1);
error = zeros(K-1,1);
z = x0;
t1 = 1;
r = 4;
s = 1 / L ;
x_k_1 = x0;
y_k_1 = x0;

for j = 1: K
    beta = (j-1)/(j+r-1);
    x_k_2 = y_k_1 - s * (A' * A * y_k_1 - A' * b + la * sign(y_k_1)) ;
    y_k_2 = x_k_2 + beta * (x_k_2 - x_k_1);
    %t(j) = j * sqrt(s);
    cost2(j) = 0.5 * norm(A * x_k_2 - b)^2 + la * norm(x_k_2, 1);
    x_k_1 = x_k_2;  
    y_k_1 = y_k_2;
end

for i = 1: K
    z = (y_k1 - L^-1 * A' * (A * y_k1 - b));
    x_k2 = max(abs(z) - la * L^-1, 0) .* sign(z);
    t2 = (1 + sqrt(1 + 4 * t1^2)) / 2;
    y_k2 = x_k1 + (t1 - 1) * t2^-1 * (x_k2 - x_k1);
    cost1(i) = 0.5 * norm(A * x_k2 - b, 2)^2 + la * norm(x_k2, 1);
    t1 = t2;
    x_k1 = x_k2;
    y_k1 = y_k2;
end

iteration = linspace(1,K,K);

figure(1)
semilogy(iteration,cost1,'color','red')
hold on
semilogy(iteration,cost2,'color','blue')
xlabel('iterations')
ylabel('f-f*')


% A' * (A * x_k1 - b)
