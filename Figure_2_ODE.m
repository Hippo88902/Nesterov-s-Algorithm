clc
clear
x0 = [1,0]';
tspan = [0.1 50];
[t,x] = ode45(@(t,x) odefcn1(t,x), tspan, x0);

error = zeros(length(x),1);
x_ = x( :, 2);

k = 1;
while (k<length(t))
    error(k) = t(k)^2 * 0.5 * (x_(k)^2); 
    k = k+1;
end
figure(1)
plot(t',error)
xlabel('t')
ylabel('t^2(f - f*)')

function dxdt = odefcn1(t,x)
  dxdt = ones(2,1);
  dxdt(1) = x(2);
  dxdt(2) = (-1/t)*dxdt(1) - x(1);
end