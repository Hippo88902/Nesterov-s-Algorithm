clear
clc
%% initial conditions
% minimizing f = 0.02*x1^2 + 0.005*x2^2, starting from x0 = (1,1) 
A = [0.02 0;0 0.005];
b = [0 0]';
x0 = [1,1]';
%% s = 0.25
K1 = 2000;
s1 = 0.25;
xk = x0;
xk_1 = x0;
yk = x0;
yk_1 = x0;
x1 = zeros(K1,1);
x2 = zeros(K1,1);
time1 = ones(K1,1);
error1 = zeros(K1,1);
k = 1;
while (k < K1)
    beta = (k-1)/(k+2);
    x1(k) = xk(1);
    x2(k) = xk(2);
    time1(k) = k * sqrt(s1);
    error1(k) = 0.02 * (xk(1)^2) + 0.005*(xk(2)^2);
    xk = yk_1 - s1 * (A * yk_1 + b);
    yk = xk + beta * (xk - xk_1);
    yk_1 = yk;
    xk_1 = xk;
    k = k + 1;
end
%% s = 0.05
K2 = 4400;
s2 = 0.05;
uk = x0;
uk_1 = x0;
vk = x0;
vk_1 = x0;
x1_ = zeros(K2,1);
x2_ = zeros(K2,1);
time2 = ones(K2,1);
error2 = zeros(K2,1);
k = 1;
while (k < K2)
    beta = (k-1)/(k+2);
    x1_(k) = uk(1);
    x2_(k) = vk(2);
    time2(k) = k * sqrt(s2);
    error2(k) = 0.02*(uk(1)^2) + 0.005*(uk(2)^2);
    uk = vk_1 - s2 .* (A * vk_1 + b);
    vk = uk + beta.*(uk - uk_1);
    uk_1 = uk;
    vk_1 = vk;
    k = k + 1;
end
%% ODE result
x0_ = [1,1,0,0]';
tspan = [0.01 1000];
[t,x] = ode45(@(t,x) odefcn(t,x), tspan, x0_);
x1 = x(:,1);
x2 = x(:,2);
error = zeros(length(x),1);
k_ = 1;
while(k_<length(x))
    error(k_) = 0.02* (x1(k_)^2) + 0.005 * (x2(k_)^2);
    k_ = k_+1 ;
end
%% x1 x2 graph
figure(1)
plot(x1, x2,'.-', 'color', 'green')
hold on
plot(x1_, x2_,'.-', 'color', 'blue')
hold on
plot(x(:,1),x(:,2),'r')
xlabel('x1')
ylabel('x2')
legend('s = 0,25','s = 0.05','ODE')
%% Error graph
figure(2)
semilogy(time1,error1,'r')
hold on
semilogy(time2,error2,'b')
hold on
semilogy(t, error,'.')
legend('s = 0,25','s = 0.05','ODE')
xlabel('t')
ylabel('f - f*')
%% ODE function
function dxdt = odefcn(t,x)
  dxdt = zeros(4,1);
  dxdt(1) = x(3);
  dxdt(2) = x(4);
  dxdt(3) = (-3/t) * dxdt(1) - 0.02 * x(1);
  dxdt(4) = (-3/t) * dxdt(2) - 0.005 * x(2);
end


