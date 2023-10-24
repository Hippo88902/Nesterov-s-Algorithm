clear
clc
x0_ = [2,0,0,0]';

tspan = [0.1 1000];

[t,x] = ode45(@(t,x) odefcn1(t,x), tspan, x0_);
error = zeros(length(x),1);
x1 = x(:,1);
x2 = x(:,2);

figure(1)
plot(x(:,1),x(:,2),'r')

figure(2)
plot(x(:,3),x(:,4),'b')


function dxdt = odefcn1(t,x)

  r = 3;
  dxdt = zeros(4,1);
  dxdt(1) = x(3);
  dxdt(2) = x(4);
  A = [20 6;6 3];
  x_g = A* [x(1) x(2)]' - [8 6]';
  
  if and(x(1) < 0, x(2) > 0)
    dxdt(3) = (-r/t)*x(3) - x_g(1) + 1;
    dxdt(4) = (-r/t)*x(4) - x_g(2) - 1;
  elseif and(x(1) > 0, x(2) < 0)
    dxdt(3) = (-r/t)*x(3) - x_g(1) - 1;
    dxdt(4) = (-r/t)*x(4) - x_g(2) + 1;
  elseif and(x(1) < 0, x(2) < 0)
    dxdt(3) = (-r/t)*x(3) - x_g(1) + 1;
    dxdt(4) = (-r/t)*x(4) - x_g(2) + 1;    
  else
    dxdt(3) = (-r/t)*x(3) - x_g(1) - 1;
    dxdt(4) = (-r/t)*x(4) - x_g(2) - 1;  
  end
  
end