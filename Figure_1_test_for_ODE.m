%A = [0.02 0;0 0.005];
%b = [0 0]';
clear
clc
x0 = [1,1,0,0]';

tspan = [1 1000];

[t,x] = ode45(@(t,x) odefcn1(t,x), tspan, x0);
error = zeros(length(x),1);
x1 = x(:,1);
x2 = x(:,2);
k = 1;

while(k<length(x))
    error(k) = 0.02* (x1(k)^2) + 0.005 * (x2(k)^2);
    k = k+1 ;
end
iteration = linspace(1, length(t), length(t));
figure(1)
plot(x(:,1),x(:,2),'r')
hold on
figure(2)
semilogy(iteration,t,'b')

function dxdt = odefcn1(t,x)
  dxdt = zeros(4,1);
  dxdt(1) = x(3);
  dxdt(2) = x(4);
  dxdt(3) = (-3/t)*dxdt(1) - 0.02.*x(1);
  dxdt(4) = (-3/t)*dxdt(2) - 0.005.*x(2);
end
