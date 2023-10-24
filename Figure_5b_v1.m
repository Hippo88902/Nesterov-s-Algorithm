clear 
clc
rng(1)
m = 500; % The number of data (row)
n = 500; % The number of variables (column)
A = normrnd(0, 1, [m,n]); % generate data
b = normrnd(0, 9, [m,1]); % generate random coefficient
la = 4;
x0 = rand(n,1);
K = 1000;
L = max(eig(A'*A));
s = 1 / L;
x_k1 = x0;
y_k1 = x0;
t = ones(K,1);
cost1 = zeros(K-1,1);
cost2 = zeros(K-1,1);
error = zeros(K-1,1);
r = 6;

for i = 1: K
    beta = (i-1)/(i+r-1);
    x_k2 = y_k1 - s * ( A' * (A * y_k1 - b) + la * sign(y_k1) ) ;
    y_k2 = x_k2 + beta * (x_k2 - x_k1);
    cost1(i) = 0.5 * norm(A * x_k2 - b)^2 + la * norm(x_k2, 1);
    x_k1 = x_k2;  
    y_k1 = y_k2;
end

x_k_1 = x0;
y_k_1 = x0;
for j = 1: K
    beta = (j-1)/(j+r-1);
    x_k_2 = y_k_1 - s *  ( A' * (A * y_k_1 - b) + la * sign(y_k_1) ) ;
    y_k_2 = x_k_2 + beta * (x_k_2 - x_k_1);
    t(j) = j * sqrt(s);
    cost2(j) = 0.5 * norm(A * x_k_2 - b)^2 + la * norm(x_k_2, 1);
    error(j) = abs(0.5 * norm(A * x_k_2 - b)^2 + la * norm(x_k_2, 1) - (0.5 * norm(A * x_k1 - b)^2 + la * norm(x_k1, 1)));
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
