clc
clear
x0_ = [1,1]';
tspan = [1 30];
[t,x] = ode45(@(t,x) odefcn1(t,x), tspan, x0_);
error_ = zeros(length(x),1);
x_ = x( :, 2);
i = linspace(1,length(x),length(x));
k_ = 1;

while (k_<length(x))
    error_(k_) = t(k_)^2 * 0.5 * (x_(k_)^2); 
    k_ = k_+1;
end

figure(1)
plot(i',error_)
xlabel('t')
ylabel('t^2(f - f*)')

function dxdt = odefcn1(t,x)
  dxdt = ones(2,1);
  dxdt(1) = x(2);
  dxdt(2) = (-2/t)*dxdt(1) - x(1);
end