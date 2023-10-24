clear
clc
a = [0 4 2 ; 1 1 1]';
A = [20 6;6 3];
b = [-8 -6]';
s1 = 0.01;
s2 = 0.001;
s3 = 0.0001;
b1 = [1 1]';
b2 = [-1 1]';
b3 = [-1 -1]';
b4 = [1 -1]';
x0 = [2,0]';
x_k1 = x0;
y_k1 = x0;
x_k_1 = x0;
y_k_1 = x0;
x_k1_ = x0;
y_k1_ = x0;
x_k2 = x0;
y_k2 = x0;
x_k_2 = x0;
y_k_2 = x0;
x_k2_ = x0;
y_k2_ = x0;

K = 10000;
xk1 = zeros(K-1,1);
xk2 = zeros(K-1,1);
xk_1 = zeros(K-1,1);
xk_2 = zeros(K-1,1);
xk1_ = zeros(K-1,1);
xk2_ = zeros(K-1,1);

k = 1;
while (k < K)
     beta = (k-1)/(k+2);
     if and(x_k2(1) > 0 , x_k2(2)>0)
        x_k2 = y_k1 - s1 * (A * y_k1 + (b + b1));
        y_k2 = x_k2 + beta * (x_k2 - x_k1);        
     
     elseif and(x_k2(1) < 0, x_k2(2) > 0)   
        x_k2 = y_k1 - s1 * (A * y_k1 + (b + b2));
        y_k2 = x_k2 + beta * (x_k2 - x_k1);

     elseif and(x_k2(1) < 0, x_k2(2) < 0)   
        x_k2 = y_k1 - s1 * (A * y_k1 + (b + b3));
        y_k2 = x_k2 + beta * (x_k2 - x_k1); 
        
     elseif and(x_k2(1) > 0, x_k2(2) < 0)
        x_k2 = y_k1 - s1 * (A * y_k1 + (b + b4));
        y_k2 = x_k2 + beta * (x_k2 - x_k1);

     elseif and(x_k2(1) == 0, x_k2(2) > 0)
        x_k2 = y_k1 - s1 * (A * y_k1 + (b + [0 1]'));
        y_k2 = x_k2 + beta * (x_k2 - x_k1);
        
     elseif and(x_k2(1) == 0, x_k2(2) < 0)
        x_k2 = y_k1 - s1 * (A * y_k1 + (b + [0 -1]'));
        y_k2 = x_k2 + beta * (x_k2 - x_k1);
        
     elseif and(x_k2(1) > 0, x_k2(2) == 0)
        x_k2 = y_k1 - s1 * (A * y_k1 + (b + [1 0]'));
        y_k2 = x_k2 + beta * (x_k2 - x_k1);
        
     elseif and(x_k2(1) < 0, x_k2(2) == 0)
        x_k2 = y_k1 - s1 * (A * y_k1 + (b + [-1 0]'));
        y_k2 = x_k2 + beta * (x_k2 - x_k1);
     
     else
        x_k2 = y_k1 - s1 * (A * y_k1 + b);
        y_k2 = x_k2 + beta * (x_k2 - x_k1);
        
     end
     
    xk1(k) = x_k2(1);
    xk2(k) = x_k2(2);
    y_k1 = y_k2;
    x_k1 = x_k2;      
    k = k + 1;
end
k1 = 1;
while (k1 < K)
     
     beta = (k1-1)/(k1+2);
     if and(x_k_2(1) > 0 , x_k_2(2)>0)
        x_k_2 = y_k_1 - s2 * (A * y_k_1 + (b + b1));
        y_k_2 = x_k_2 + beta * (x_k_2 - x_k_1); 

     elseif and(x_k_2(1) < 0, x_k_2(2) > 0)   
        x_k_2 = y_k_1 - s2 * (A * y_k_1 + (b + b2));
        y_k_2 = x_k_2 + beta * (x_k_2 - x_k_1); 

     elseif and(x_k_2(1) < 0, x_k_2(2) < 0)   
        x_k_2 = y_k_1 - s2 * (A * y_k_1 + (b + b3));
        y_k_2 = x_k_2 + beta * (x_k_2 - x_k_1); 

     elseif and(x_k_2(1) > 0, x_k_2(2) < 0)
        x_k_2 = y_k_1 - s2 * (A * y_k_1 + (b + b4));
        y_k_2 = x_k_2 + beta * (x_k_2 - x_k_1); 

     elseif and(x_k_2(1) == 0, x_k_2(2) > 0)
        x_k_2 = y_k_1 - s2 * (A * y_k_1 + (b + [0 1]'));
        y_k_2 = x_k_2 + beta * (x_k_2 - x_k_1); 

     elseif and(x_k_2(1) == 0, x_k_2(2) < 0)
        x_k_2 = y_k_1 - s2 * (A * y_k_1 + (b + [0 -1]'));
        y_k_2 = x_k_2 + beta * (x_k_2 - x_k_1); 
        
     elseif and(x_k_2(1) > 0, x_k_2(2) == 0)
        x_k_2 = y_k_1 - s2 * (A * y_k_1 + (b + [1 0]'));
        y_k_2 = x_k_2 + beta * (x_k_2 - x_k_1);   
        
     elseif and(x_k_2(1) < 0, x_k_2(2) == 0)
        x_k_2 = y_k_1 - s2 * (A * y_k_1 + (b + [-1 0]'));
        y_k_2 = x_k_2 + beta * (x_k_2 - x_k_1);  
     
     else
        x_k_2 = y_k_1 - s2 * (A * y_k_1 + b);
        y_k_2 = x_k_2 + beta * (x_k_2 - x_k_1); 
        
     end
    
    xk_1(k1) = x_k_2(1);
    xk_2(k1) = x_k_2(2);
    y_k_1 = y_k_2;
    x_k_1 = x_k_2;     
    k1 = k1 + 1;

end
k2 = 1;
while (k2 < K)
     
     beta = (k2-1)/(k2+2);
     if and(x_k2_(1) > 0 , x_k2_(2)>0)
        x_k2_ = y_k1_ - s3 * (A * y_k1_ + (b + b1));
        y_k2_ = x_k2_ + beta * (x_k2_ - x_k1_);        
     
     elseif and(x_k2_(1) < 0, x_k2_(2) > 0)    
        x_k2_ = y_k1_ - s3 * (A * y_k1_ + (b + b2));
        y_k2_ = x_k2_ + beta * (x_k2_ - x_k1_); 
        
     elseif and(x_k2_(1) < 0, x_k2_(2) < 0)    
        x_k2_ = y_k1_ - s3 * (A * y_k1_ + (b + b3));
        y_k2_ = x_k2_ + beta * (x_k2_ - x_k1_); 
        
     elseif and(x_k2_(1) > 0, x_k2_(2) < 0) 
        x_k2_ = y_k1_ - s3 * (A * y_k1_ + (b + b4));
        y_k2_ = x_k2_ + beta * (x_k2_ - x_k1_); 
        
     elseif and(x_k2_(1) == 0, x_k2_(2) > 0) 
        x_k2_ = y_k1_ - s3 * (A * y_k1_ + (b + [0 1]'));
        y_k2_ = x_k2_ + beta * (x_k2_ - x_k1_); 
        
     elseif and(x_k2_(1) == 0, x_k2_(2) < 0) 
        x_k2_ = y_k1_ - s3 * (A * y_k1_ + (b + [0 -1]'));
        y_k2_ = x_k2_ + beta * (x_k2_ - x_k1_); 
        
     elseif and(x_k2_(1) > 0, x_k2_(2) == 0) 
        x_k2_ = y_k1_ - s3 * (A * y_k1_ + (b + [1 0]'));
        y_k2_ = x_k2_ + beta * (x_k2_ - x_k1_);  
        
     elseif and(x_k2_(1) < 0, x_k2_(2) == 0) 
        x_k2_ = y_k1_ - s3 * (A * y_k1_ + (b + [-1 0]'));
        y_k2_ = x_k2_ + beta * (x_k2_ - x_k1_); 
     
     else 
        x_k2_ = y_k1_ - s3 * (A * y_k1_ + b);
        y_k2_ = x_k2_ + beta * (x_k2_ - x_k1_);  
     end
    
    xk1_(k2) = x_k2_(1);
    xk2_(k2) = x_k2_(2); 
    y_k1_ = y_k2_;
    x_k1_ = x_k2_;         
    k2 = k2 + 1;

end

% test for numerical convergence 
test = a' * ( a*y_k1 - [4 2 0]' ) + [-1 1]';

x0_ = [2,0,0,0]';
tspan = [0.1 10];
options = odeset('RelTol',1e-8,'Stats','on','OutputFcn',@odeplot);
[t,x] = ode45(@(t,x) odefcn1(t,x), tspan, x0_,options);
x1 = x(:,1);
x2 = x(:,2);
% eigen vector of matrix A
[V, D] = eig(A);
v1 = D(1,1)*V(:,1);
v2 = D(2,2)*V(:,2);
% The comparision of different s and ODE 
figure(1)
% s = 0.01
plot(xk1,xk2,'-.','color','red')
hold on
plot(xk1(ceil(1/sqrt(s1))),xk2(ceil(1/sqrt(s1))),'ro')
hold on
plot(xk1(ceil(2/sqrt(s1))),xk2(ceil(2/sqrt(s1))),'r+')
hold on
plot(xk1(K-1),xk2(K-1),'r^')
hold on
% s = 0.001
plot(xk_1,xk_2,'-.','color','blue')
hold on
plot(xk_1(ceil(1/sqrt(s2))),xk_2(ceil(1/sqrt(s2))),'bo')
hold on
plot(xk_1(ceil(2/sqrt(s2))),xk_2(ceil(2/sqrt(s2))),'b+')
hold on
plot(xk_1(K-1),xk_2(K-1),'b^')
hold on
% s = 0.0001
plot(xk1_,xk2_,'-.','color','black')
hold on
plot(xk1_(ceil(1/sqrt(s3))),xk2_(ceil(1/sqrt(s3))),'o','color','black')
hold on
plot(xk1_(ceil(2/sqrt(s3))),xk2_(ceil(2/sqrt(s3))),'+','color','black')
hold on
plot(xk1_(K-1),xk2_(K-1),'^','color','black')
hold on
%% ODE
plot(x(:,1),x(:,2),'yellow')
hold on
quiver(xk1(K-1), xk2(K-1), v1(1), v1(2),'g')
xlabel('x1')
ylabel('x2')
error1 = zeros(K-1,1);
error2 = zeros(K-1,1);
error3 = zeros(K-1,1);

for i = 2:K-1
    error1(i) = xk1(i) - xk1(i-1);
    error2(i) = xk2(i) - xk2(i-1);
end 

y = [4 2 0]';
x_input = zeros(K-1,2);
x_input(:,1) = xk1;
x_input(:,2) = xk2;
for i = 1:K-1
    error3(i) = norm(y-a*x_input(i,:)',2)^2/2 + norm(x_input(i,:)',1) - norm(y)^2/2 ;
end


iteration = linspace(1,500,500);
figure(2)
semilogy(iteration,abs(error1(1:500)))
xlabel('iteration')
ylabel('x1 error')
figure(3)
semilogy(iteration,abs(error2(1:500)))
xlabel('iteration')
ylabel('x2 error')
%figure(4)
%semilogy(iteration,abs(error3(1:500)))


function dxdt = odefcn1(t,x)
  r = 3;
  dxdt = zeros(4,1);
  dxdt(1) = x(3);
  dxdt(2) = x(4);
  A = [20 6;6 3];
  x_ = A* [x(1) x(2)]' - [8 6]';
  if and(x(1) < 0, x(2) > 0)
    dxdt(3) = (-r/t)*dxdt(1) - x_(1) + 1;
    dxdt(4) = (-r/t)*dxdt(2) - x_(2) - 1;
  
  elseif and(x(1) > 0, x(2) < 0)
    dxdt(3) = (-r/t)*dxdt(1) - x_(1) - 1;
    dxdt(4) = (-r/t)*dxdt(2) - x_(2) + 1;
  
  elseif and(x(1) < 0, x(2) < 0)
    dxdt(3) = (-r/t)*dxdt(1) - x_(1) + 1;
    dxdt(4) = (-r/t)*dxdt(2) - x_(2) + 1;    
  
  else
    dxdt(3) = (-r/t)*dxdt(1) - x_(1) - 1;
    dxdt(4) = (-r/t)*dxdt(2) - x_(2) - 1;  
  end
  
end

