clear
clc
x0 = [-0.1,0]';
K = 501;
s = 0.05;
x_k1 = x0;
y_k1 = x0;
x_k_1 = x0;
x_k2 = x0;
x_k_2 = x0;
y_kg = [1 1]';
x_kg = [1 1]'; 
xk1 = zeros(K-1,1);
xk2 = zeros(K-1,1);
gradient1 = zeros(K-1,2);
k = 1;
r = 5;
% Nesterov's scheme
while (k < K)
    
    if (y_k1(2) == 0)
        y_kg(1) = 3 * y_k1(1)^3/abs(y_k1(1)) + 0.002 * y_k1(1);
        y_kg(2) = 0.002 * y_k1(1);  
    elseif (y_k1(1) == 0)
        y_kg(1) = 0.002 * y_k1(2);
        y_kg(2) = 15 * y_k1(2)^3/abs(y_k1(2)) + 0.002 * y_k1(2);            
    else
        y_kg(1) = 3 * y_k1(1)^3/abs(y_k1(1)) + 0.002 * (y_k1(1) + y_k1(2));
        y_kg(2) = 15 * y_k1(2)^3/abs(y_k1(2)) + 0.002 * (y_k1(1) + y_k1(2));        
    end    
    
    beta = (k-1)/(k+r-1);
    x_k2 = y_k1 - s * y_kg ;
    y_k2 = x_k2 + beta * (x_k2 - x_k1);
  
    xk1(k) = x_k2(1); %紀錄每次迭代的x1值
    xk2(k) = x_k2(2); %紀錄每次迭代的x2值
    gradient1(k,:) = y_kg; %紀錄每次迭代的gradient_f值
    
    y_k1 = y_k2;
    x_k1 = x_k2;
    k = k + 1;
end
j = 1;
J = 100000;
xk_1 = zeros(J-1,1);
xk_2 = zeros(J-1,1);
gradient2 = zeros(J-1,2);
% Gradient descent
while (j<J)
    
    if (x_k_1(2) == 0)
        x_kg(1) = 3 * x_k_1(1)^3/abs(x_k_1(1)) + 0.002 * x_k_1(1);
        x_kg(2) = 0.002 * x_k_1(1);  
    elseif (y_k1(1) == 0)
        x_kg(1) = 0.002 * x_k_1(2);
        x_kg(2) = 15 * x_k_1(2)^3/abs(x_k_1(2)) + 0.002 * x_k_1(2);
    else    
        x_kg(1) = 3 * x_k_1(1)^3/abs(x_k_1(1)) + 0.002 * (x_k_1(1) + x_k_1(2));
        x_kg(2) = 15 * x_k_1(2)^3/abs(x_k_1(2)) + 0.002 * (x_k_1(1) + x_k_1(2));    
    end
    
    x_k_2 = x_k_1 - s * x_kg ;
    
    xk_1(j) = x_k_2(1); %紀錄每次迭代的x1值
    xk_2(j) = x_k_2(2); %紀錄每次迭代的x2值
    gradient2(j,:) = x_kg; %紀錄每次迭代的gradient_f值
    
    x_k_1 = x_k_2;
    j = j + 1;
    
end

tspan1 = [0.1 80];
x0_ = [-0.1,0,0,0]';
[t,x] = ode45(@(t,x) odefcn1(t,x), tspan1, x0_);
x1 = x(:,1);
x2 = x(:,2);

tspan2 = [0.1 1000];
x_0 = [-0.1,0]';
[t2,X] = ode45(@(t2,X) odefcn2(X), tspan2, x_0);
X1 = X(:,1);
X2 = X(:,2);


figure(1)
plot(x(:,1),x(:,2),'g')
hold on
plot(X(:,1),X(:,2),'y')
hold on
plot(xk1,xk2,'r')
hold on
plot(xk_1,xk_2,'b')
xlabel('x1')
ylabel('x2')


function dxdt = odefcn1(t,x)
  r = 3;
  dxdt = zeros(4,1);
  dxdt(1) = x(3);
  dxdt(2) = x(4);
  dxdt(3) = (-r/t)*x(3) - (3 * sign(x(1)) * x(1)^2 + 0.002 * (x(1) + x(2)));
  dxdt(4) = (-r/t)*x(4) - (15 * sign(x(2)) * x(2)^2 + 0.002 * (x(1) + x(2)));
end

function dxdt = odefcn2(x)
  dxdt = zeros(2,1);
  dxdt(1) = - (3 * sign(x(1)) * x(1)^2 + 0.002 * (x(1) + x(2)));
  dxdt(2) = - (15 * sign(x(2)) * x(2)^2 + 0.002 * (x(1) + x(2)));
end
