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
    
    beta = (k-1)/(k+2);
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

figure(1)
plot(xk1,xk2,'r')
hold on
plot(xk1(1),xk2(1),'ro')
hold on
plot(xk1(10),xk2(10),'ro')
hold on
plot(xk1(20),xk2(20),'ro')
hold on
plot(xk1(30),xk2(30),'ro')
hold on
plot(xk1(45),xk2(45),'ro')
hold on
plot(xk1(60),xk2(60),'ro')
hold on
plot(xk1(90),xk2(90),'ro')
hold on
plot(xk1(120),xk2(120),'ro')
hold on
plot(xk1(150),xk2(150),'ro')
hold on
plot(xk1(190),xk2(190),'ro')
hold on
plot(xk1(250),xk2(250),'ro')
hold on
plot(xk1(300),xk2(300),'ro')
hold on
plot(xk_1,xk_2,'b')
hold on
plot(xk_1(1),xk_2(1),'bo')
hold on
plot(xk_1(10),xk_2(10),'bo')
hold on
plot(xk_1(20),xk_2(20),'bo')
hold on
plot(xk_1(30),xk_2(30),'bo')
hold on
plot(xk_1(45),xk_2(45),'bo')
hold on
plot(xk_1(60),xk_2(60),'bo')
hold on
plot(xk_1(90),xk_2(90),'bo')
hold on
plot(xk_1(120),xk_2(120),'bo')
hold on
plot(xk_1(150),xk_2(150),'bo')
hold on
plot(xk_1(190),xk_2(190),'bo')
hold on
plot(xk_1(250),xk_2(250),'bo')
hold on
plot(xk_1(300),xk_2(300),'bo')
hold on
xlabel('x1')
ylabel('x2')
