clear
clc
x0 = 1;
K = 5000;
s = 0.0001;
xk = x0;
xk_1 = x0;
yk = x0;
yk_1 = x0;
x = zeros(K-1,1);
t = ones(K-1,1);
error = zeros(K-1,1);
r=1;
k = 1;
while (k < K)
    beta = (k-1)/(k+r-1);
    xk = yk_1 - s * yk_1;
    yk = xk + beta*(xk - xk_1);
    x(k) = xk;
    t(k) = k * sqrt(s);
    error(k) = t(k)^2 * (0.5 * xk^2);
    yk_1 = yk;
    xk_1 = xk;     
    k = k + 1;
end

iteration = linspace(1,K-1,K-1);
figure(1)
plot(iteration, error, '-.', 'color', 'red')
xlabel('iteration')
ylabel('sk^2(f - f*)')


