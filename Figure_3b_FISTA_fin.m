clear
clc
A = [1 1];
K = 1000;
x0 = [-0.1, 0]';
xk = x0;
xk_1 = x0;
yk = x0;
L = max(eig(A' * A)) * 0.493;
t1 = 1;
x1 = zeros(K,1);
x2 = zeros(K,1);
error = zeros(K-1,1);

for i = 1: K
    g = abs(xk(1))^3 + 5 * abs(xk(2))^3;
    f = 0.001 * (xk(1) + xk(2))^2;
    x1(i) = xk(1);
    x2(i) = xk(2);
    error(i) = f + g;
    c1 = L * yk(1) -  0.002 * (yk(1) + yk(2));
    c2 = L * yk(2) -  0.002 * (yk(1) + yk(2));    
    % x1 
    if (c1 > 0)
        xk(1) =  (-L + sqrt(L^2 + 12 * L * c1)) / 6 ;
    elseif (c1 < 0)
        xk(1) = (L - sqrt(L^2 - 12 * L * c1)) / 6 ;
    else
        xk(1) = 0;
    end
    % x2
    if (c2 > 0)
        xk(2) = (-L + sqrt(L^2 + 60 * L * c2)) / 30 ;
    elseif (c2 < 0)
        xk(2) = (L - sqrt(L^2 - 60 * L * c2)) / 30 ;
    else
        xk(2) = 0;
    end
    t2 = (1 + sqrt(1 + 4 * t1^2)) / 2;
    yk1 = xk + (t1 - 1) * t2^-1 * (xk - xk_1);
    t1 = t2;
    xk_1 = xk;
    yk = yk1;
end 

figure(1)
plot(x1,x2,'color','blue')
xlabel('x1')
ylabel('x2')

%figure(2)
%iteration = linspace(1,K,K);
%semilogy(iteration, error, 'color','red')
%xlabel('iterations')
%ylabel('error')