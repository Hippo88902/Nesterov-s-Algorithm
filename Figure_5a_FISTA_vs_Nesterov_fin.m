clear clc
rng(10)
m = 100; % The number of data (row)
n = 500; % The number of variables (column)
A = normrnd(0, 1, [m,n]); % generate data
b = normrnd(0, 25, [m,1]); % generate random coefficient
la = 4;
x0 = 0.1 * ones(n,1);
K = 5000;
L = max(eig(A' * A));
s = 1 / L;
t = ones(K,1);
cost = zeros(K-1,1);
error1 = zeros(K-1,1);
error2 = zeros(K-1,1);
error3 = zeros(K-1,1);
% FISTA
xk = x0;
xk_1 = x0;
yk = x0;
yk1 = x0;
z = x0;
t1 = 1;
for i = 1: 30000
    z = (yk1 - L^-1 * A' * (A * yk1 - b));
    xk = max(abs(z) - la * L^-1, 0) .* sign(z);
    t2 = (1 + sqrt(1 + 4 * t1^2)) / 2;
    yk1 = xk + (t1 - 1) * t2^-1 * (xk - xk_1);
    cost(i) = 0.5 * norm(A * xk - b, 2)^2 + la * norm(xk, 1);
    t1 = t2;
    xk_1 = xk;
    yk = yk1;
end
f_opt = 0.5 * norm(A * xk - b)^2 + la * norm(xk, 1);
% Nesterov
uk = x0;
uk_1 = x0;
vk = x0;
vk_1 = x0;
for i = 1: K
    beta = (i-1)/(i+2);
    u = vk_1 - s * A' *(A * vk_1 - b);
    argmin = max( abs(u) - la * s, 0 ) .* sign(u);
    Gs = (vk_1 - argmin) / s;
    uk = vk_1 - s * Gs ;
    vk = uk + beta * (uk - uk_1);
    error1(i) = abs(0.5 * norm(A * uk - b)^2 + la * norm(uk, 1) - f_opt);
    uk_1 = uk;  
    vk_1 = vk;
end

mk = x0;
mk_1 = x0;
nk = x0;
nk_1 = x0;
for i = 1: K
    beta = (i-1)/(i+3);
    u = nk_1 - s * A' *(A * nk_1 - b);
    argmin = max( abs(u) - la * s, 0 ) .* sign(u);
    Gs = (nk_1 - argmin) / s;
    mk = nk_1 - s * Gs ;
    nk = mk + beta * (mk - mk_1);
    error2(i) = abs(0.5 * norm(A * mk - b)^2 + la * norm(mk, 1) - f_opt);
    mk_1 = mk;  
    nk_1 = nk;
end

sk = x0;
sk_1 = x0;
tk = x0;
tk_1 = x0;
for i = 1: K
    beta = (i-1)/(i+4);
    u = tk_1 - s * A' *(A * tk_1 - b);
    argmin = max( abs(u) - la * s, 0 ) .* sign(u);
    Gs = (tk_1 - argmin) / s;
    sk = tk_1 - s * Gs ;
    tk = sk + beta * (sk - sk_1);
    error3(i) = abs(0.5 * norm(A * sk - b)^2 + la * norm(sk, 1) - f_opt);
    sk_1 = sk;  
    tk_1 = tk;
end


iteration1 = linspace(1,30000,30000);
figure(1)
plot(iteration1,cost,'color','red')
xlabel('iterations')
ylabel('cost')

iteration2 = linspace(1,K,K);
figure(2)
semilogy(iteration2,error1,'color','black')
hold on
semilogy(iteration2,error2,'color','red')
hold on
semilogy(iteration2,error3,'color','blue')
xlabel('iterations')
ylabel('f-f*')
legend('r = 3','r = 4','r = 5')

% A' * (A * xk - b)
