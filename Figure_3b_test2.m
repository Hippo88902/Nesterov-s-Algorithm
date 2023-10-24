clear
clc
x0 = [-0.01,0]';
K = 300;
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
while (k < K)
    
    if and(y_k1(1) < 0, y_k1(2) < 0)
        y_kg(1) = -3 * y_k1(1)^2 + 0.002 * (y_k1(1) + y_k1(2));
        y_kg(2) = -15 * y_k1(2)^2 + 0.002 * (y_k1(1) + y_k1(2));
        
    elseif and(y_k1(1) < 0, y_k1(2) > 0)
        y_kg(1) = -3 * y_k1(1)^2 + 0.002 * (y_k1(1) + y_k1(2));
        y_kg(2) = 15 * y_k1(2)^2 + 0.002 * (y_k1(1) + y_k1(2));    
       
    elseif and(y_k1(1) > 0, y_k1(2) < 0)
        y_kg(1) = 3 * y_k1(1)^2 + 0.002 * (y_k1(1) + y_k1(2));
        y_kg(2) = -15 * y_k1(2)^2 + 0.002 * (y_k1(1) + y_k1(2));
  
    elseif and(y_k1(1) > 0, y_k1(2) > 0)
        y_kg(1) = 3 * y_k1(1)^2 + 0.002 * (y_k1(1) + y_k1(2));
        y_kg(2) = 15 * y_k1(2)^2 + 0.002 * (y_k1(1) + y_k1(2));    

    elseif and(y_k1(1) < 0, y_k1(2) == 0)
        y_kg(1) = -3 * y_k1(1)^2 + 0.002 * y_k1(1);
        y_kg(2) = 0.002 * y_k1(1);  
        
    elseif and(y_k1(1) > 0, y_k1(2) == 0)
        y_kg(1) = 3 * y_k1(1)^2 + 0.002 * y_k1(1);
        y_kg(2) = 0.002 * y_k1(1);
    
    elseif and(y_k1(1) == 0, y_k1(2) < 0)
        y_kg(1) = 0.002 * y_k1(2);
        y_kg(2) = -15 * y_k1(2)^2 + 0.002 * y_k1(2);        

    elseif and(y_k1(1) == 0, y_k1(2) > 0)
        y_kg(1) = 0.002 * y_k1(2);
        y_kg(2) = 15 * y_k1(2)^2 + 0.002 * y_k1(2);        
        
    end
    
    beta = (k-1)/(k+2);
    x_k2 = y_k1 - s * y_kg ;
    y_k2 = x_k2 + beta * (x_k2 - x_k1);
    xk1(k) = x_k2(1);
    xk2(k) = x_k2(2);
    gradient1(k,:) = y_kg;
    y_k1 = y_k2;
    x_k1 = x_k2;
    k = k + 1;

end

j = 1;
J = 10000;
xk_1 = zeros(J-1,1);
xk_2 = zeros(J-1,1);
gradient2 = zeros(J-1,2);
while (j<J)
    
    if and(x_k_1(1) < 0, x_k_1(2) < 0)
        x_kg(1) = -3 * x_k_1(1)^2 + 0.002 * (x_k_1(1) + x_k_1(2));
        x_kg(2) = -15 * x_k_1(2)^2 + 0.002 * (x_k_1(1) + x_k_1(2));
    
    elseif and(x_k_1(1) < 0, x_k_1(2) > 0)    
        x_kg(1) = -3 * x_k_1(1)^2 + 0.002 * (x_k_1(1) + x_k_1(2));
        x_kg(2) = 15 * x_k_1(2)^2 + 0.002 * (x_k_1(1) + x_k_1(2));
    
    elseif and(x_k_1(1) > 0, x_k_1(2) < 0)
        x_kg(1) = 3 * x_k_1(1)^2 + 0.002 * (x_k_1(1) + x_k_1(2));
        x_kg(2) = -15 * x_k_1(2)^2 + 0.002 * (x_k_1(1) + x_k_1(2));
    
    else
        x_kg(1) = 3 * x_k_1(1)^2 + 0.002 * (x_k_1(1) + x_k_1(2));
        x_kg(2) = 15 * x_k_1(2)^2 + 0.002 * (x_k_1(1) + x_k_1(2));
   
    end
    
    x_k_2 = x_k_1 - s * x_kg ;
    xk_1(j) = x_k_2(1);
    xk_2(j) = x_k_2(2);
    gradient2(j,:) = x_kg;
    x_k_1 = x_k_2;
    j = j + 1;
    
end

figure(1)
plot(xk1,xk2,'r')
hold on
plot(xk1(5),xk2(5),'o')
hold on
plot(xk_1,xk_2,'b')
xlabel('x1')
ylabel('x2')
