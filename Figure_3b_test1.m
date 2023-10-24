tspan2 = [0.1 1000];
x_0 = [-0.1,0]';
[t2,X] = ode45(@(t2,X) odefcn2(X), tspan2, x_0);
X1 = X(:,1);
X2 = X(:,2);

figure(1)
plot(X(:,1),X(:,2),'r')

function dxdt = odefcn2(x)
  dxdt = zeros(2,1);
  dxdt(1) = - (3 * sign(x(1)) * x(1)^2 + 0.002 * (x(1) + x(2)));
  dxdt(2) = - (15 * sign(x(2)) * x(2)^2 + 0.002 * (x(1) + x(2)));
end