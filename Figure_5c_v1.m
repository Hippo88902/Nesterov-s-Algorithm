clear 
clc
rng(1)
m = 100; % The number of data (row)
n = 500; % The number of variables (column)
A = normrnd(0, 1, [m,n]); % generate data
b = normrnd(0, 25, [m,1]); % generate random coefficient
la = 4;
x0 = 0.1 * ones(n,1);
K = 10;
L = max(eig(A'*A));
s = 1 / L;
x_k1 = x0;
y_k1 = x0;
t = ones(K,1);
cost1 = zeros(K-1,1);
cost2 = zeros(K-1,1);
error = zeros(K-1,1);
r = 3;

for i = 1: K
    beta = (i-1)/(i+r-1);
    x_k2 = y_k1 - 2 * s *  (A' * (A * y_k1 - b));
    y_k2 = x_k2 + beta * (x_k2 - x_k1);
    cost1(i) = norm(A * x_k2 - b)^2;
    x_k1 = x_k2;  
    y_k1 = y_k2;
end

x_k_1 = x0;
y_k_1 = x0; 
for j = 1:K 
    beta = (j - 1)/(j + r - 1);
    x_k_2 = y_k_1 - 2 * s *  A' * (A * y_k_1 - b);
    y_k_2 = x_k_2 + beta * (x_k_2 - x_k_1);
    t(j) = j * sqrt(s);
    cost2(j) = norm(A * x_k_2 - b)^2;
    error(j) = abs(norm(A * x_k_2 - b)^2 - norm(A * x_k1 - b)^2);
    x_k_1 = x_k_2;  
    y_k_1 = y_k_2;
end


iteration = linspace(1,K,K);

figure(1)
semilogy(iteration,cost1,'color','red')
hold on
semilogy(iteration,cost2,'color','blue')
xlabel('iterations')
ylabel('cost')

figure(2)
semilogy(iteration,error,'color','red')
xlabel('iterations')
ylabel('f-f*')


% A' * (A * x_k1 - b)
