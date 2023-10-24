clear 
clc
A = [0 4 2 ; 1 1 1]';
la = 1;
K = 800000;
x0 = [2 0]';
b = [4 2 0]';
xk = x0;
xk_1 = x0;
yk = x0;
yk_1 = x0;
error = zeros(K-1,1);
%L = max(eig(A' * A));
L = 3000;
z = x0;
t1 = 1;
x1 = zeros(K, 1);
x2 = zeros(K, 1);

for i = 1: K
    error(i) = 0.5 * norm(b - A * xk, 2)^2 + la * norm(xk, 1);
    x1(i) = xk(1);
    x2(i) = xk(2);
    z = (yk_1 + L^-1 * A' * (b - A * yk_1));%可以換S
    xk = max(abs(z) - la * L^-1, 0) .* sign(z);
    t2 = (1 + sqrt(1 + 4 * t1^2)) / 2;
    yk = xk_1 + (t1 - 1) * t2^-1 * (xk - xk_1);    
    xk_1 = xk;
    yk_1 = yk;
    t1 = t2;
end

figure(1)
plot(x1,x2,'color','blue')
xlabel('x1')
ylabel('x2')

iteration = linspace(1,K,K);
figure(2)
semilogy(iteration,error,'color','red')
xlabel('iterations')
ylabel('f-f*')