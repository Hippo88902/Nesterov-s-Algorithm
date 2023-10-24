clear 
clc
A = [0.02 0;0 0.005];
x0 = [1,1]';
b = [0 0]';
K = 1000;
xk = x0;
xk_1 = x0;
yk = x0;
error = zeros(K-1,1);
L =  max(eig(A));
t1 = 1;
x1 = zeros(K, 1);
x2 = zeros(K, 1);
time = zeros(K, 1);
s = 0.005;

for i = 1: K
    xk = yk - s * L^-1 * (A * yk + b) ;
    t2 = (1 + sqrt(1 + 4 * t1^2)) / 2 ;
    yk1 = xk + (t1 - 1) * t2^-1 * (xk - xk_1);
    x1(i) = xk_1(1);
    x2(i) = xk_1(2);
    error(i) = 0.02 * (xk_1(1)^2) + 0.005 * (xk_1(2)^2);
    time(i) = t1;    
    xk_1 = xk;
    yk = yk1;
    t1 = t2;
end

figure(1)
plot(x1,x2,'color','blue')
xlabel('x1')
ylabel('x2')

iteration = linspace(1,K,K);
figure(2)
plot(iteration,x1,'color','red')
axis([0 K -0.01 0.01])
xlabel('iterations')
ylabel('x1')