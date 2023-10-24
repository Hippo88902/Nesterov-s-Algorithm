clear 
clc
rng(1)
n = 500; % The number of variables (column)
tspan1 = [0.1 80];
X = ones(500,1);
x0 = zeros(500,1);
x0_ = [X x0];
[t,x] = ode45(@(t,x) odefcn1(t,x), tspan1, x0_);

%x1 = x(:,1);
%x2 = x(:,2);

%figure(1)
%plot(x(:,1),x(:,2),'g')
%xlabel('x1')
%ylabel('x2')


function dxdt = odefcn1(t,x)
  m = 100; % The number of data (row)
  n = 500; % The number of variables (column)
  A = normrnd(0, 1, [m,n]); % generate data
  b = normrnd(0, 25, [m,1]); % generate random coefficient
  la = 4;
  r = 1000;
  dxdt = zeros(1000,1);
  
  for i = 1:500
      dxdt(i) = x(500 + i);
  end    
  
  gradient = A' * A * x(1:500) - A' * b + la * sign(x(1:500));
  
  for j = 1:500
      dxdt(500 + j) = (-r/t) * x(j) - gradient(j);
  end
  
end