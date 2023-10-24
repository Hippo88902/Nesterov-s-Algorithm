clear 
clc
A = [0.02 0;0 0.005];
x0 = [1,1]';
b = [0 0]';
%% Nesterov
K = 10000;
xk = x0;
xk_1 = x0;
yk = x0;
yk_1 = x0;
theta = 1;
error1 = zeros(K-1,1);
x1 = zeros(K, 1);
x2 = zeros(K, 1);
th = zeros(K, 1);
s = 0.005;
for i = 1: K
    xk = yk_1 - s * (A * yk_1 + b) ;
    theta1 = (sqrt(theta^4 + 4*theta^2) - theta^2) /2 ;
    yk = xk + (theta^-1 - 1) * theta1 * (xk - xk_1);
    x1(i) = xk_1(1);
    x2(i) = xk_1(2);
    error1(i) = 0.02 * (xk_1(1)^2) + 0.005 * (xk_1(2)^2);
    th(i) = theta;    
    theta = theta1;
    xk_1 = xk;
    yk_1 = yk;
end
%% FISTA
uk = x0;
uk_1 = x0;
vk = x0;
L =  max(eig(A));
t1 = 1;
error2 = zeros(K-1,1);
u1 = zeros(K, 1);
u2 = zeros(K, 1);
t = zeros(K, 1);
for i = 1: K
    uk = vk - 0.0001 * L^-1 * (A * vk + b) ;
    t2 = (1 + sqrt(1 + 4 * t1^2)) / 2 ;
    vk1 = uk + (t1 - 1) * t2^-1 * (uk - uk_1);
    u1(i) = uk_1(1);
    u2(i) = uk_1(2);
    error2(i) = 0.02 * (uk_1(1)^2) + 0.005 * (uk_1(2)^2);
    t(i) = t1;    
    uk_1 = uk;
    vk = vk1;
    t1 = t2;
end
%% New method
sk = x0;
sk_1 = x0;
tk = x0;
tk_1 = x0;
la = 1;
la_1 = 1;
error3 = zeros(K-1,1);
s1 = zeros(K, 1);
s2 = zeros(K, 1);
lamda = zeros(K, 1);
s = 0.005;
r = 5;
for i = 1: K
    sk = tk_1 - s * (A * tk_1 + b) ;
    la = (r-1) / (i+r-1);
    tk = sk + (la_1^-1 - 1) * la * (sk - sk_1);
    s1(i) = sk_1(1);
    s2(i) = sk_1(2);
    error1(i) = 0.02 * (sk_1(1)^2) + 0.005 * (sk_1(2)^2);
    lamda(i) = la_1;    
    la_1 = la;
    sk_1 = sk;
    tk_1 = tk;
end
%% graph
figure(1)
plot(x1,x2,'color','blue')
hold on 
plot(u1,u2,'color','red')
hold on
plot(s1,s2,'color','black')
xlabel('x1')
ylabel('x2')
legend('Nesterov','FISTA','Gs')
%% graph2
iteration = linspace(1,K,K);
figure(2)
plot(iteration,x1,'color','red')
hold on
plot(iteration,u1,'color','blue')
%axis([0 K -0.01 0.01])
xlabel('iterations')
ylabel('x1')