clear
clc
A = 0.001 * [1,1];
K = 1500;
x0 = [-0.1,0]';
x_k1 = x0;
y_k1 = x0;
x_k2 = x0;
error = zeros(K-1,1);
L = max(eig(A' * A));
t1 = 1 ;
x1 = zeros(K,1);
x2 = zeros(K,1);
for i = 1: K
    g = abs(x_k2(1))^3 + 5 * abs(x_k2(2))^3;
    f = 0.001 * (x_k2(1) + x_k2(2))^2;
    x1(i) = x_k2(1);
    x2(i) = x_k2(2);
    error(i) = f + g;
    v = 0.002 * (y_k1(1) + y_k1(2));
    gradient_f = 0.002 * [v v];
    fun = @(x_k2)abs(x_k2(1))^3 + 5 * abs(x_k2(2))^3 + L/2 * norm(x_k2 - (y_k1 - L^-1 * gradient_f))^2;
    x_k2 = fminsearch(fun,y_k1);
    t2 = (1 + sqrt(1 + 4 * t1^2)) / 2;
    y_k2 = x_k1 + (t1 - 1) * t2^-1 * (x_k2 - x_k1);
    t1 = t2;
    x_k1 = x_k2;
    y_k1 = y_k2;
end 

iteration = linspace(1,K,K);
figure(1)
semilogy(iteration,error,'color','red')
xlabel('iterations')
ylabel('error')

figure(2)
plot(x1,x2,'color','blue')
xlabel('x1')
ylabel('x2')
    
