clear
clc
a = [0 4 2 ; 1 1 1]';
A = [20 6;6 3];
b = [-8 -6]';
s1 = 0.01;
s2 = 0.001;
s3 = 0.0001;
x0 = [2,0]';
K = 10000;

k = 1;
xk = x0;
xk_1 = x0;
yk = x0;
yk_1 = x0;
x1 = zeros(K-1,1);
x2 = zeros(K-1,1);
while (k < K)
    x1(k) = xk(1);
    x2(k) = xk(2);
    beta = (k-1)/(k+2);
    xk = yk_1 - s1 * (A * yk_1 + (b + sign(xk)));
    yk = xk + beta * (xk - xk_1);            
    yk_1 = yk;
    xk_1 = xk;      
    k = k + 1;
end

k1 = 1;
uk = x0;
uk_1 = x0;
vk = x0;
vk_1 = x0;
u1 = zeros(K-1,1);
u2 = zeros(K-1,1);
while (k1 < K)    
    u1(k1) = uk(1);
    u2(k1) = uk(2);
    beta = (k1-1)/(k1+2);
    uk = vk_1 - s2 * (A * vk_1 + (b + sign(uk)));
    vk = uk + beta * (uk - uk_1); 
    uk_1 = uk;
    vk_1 = vk;     
    k1 = k1 + 1;
end

k2 = 1;
mk = x0;
mk_1 = x0;
nk = x0;
nk_1 = x0;
m1 = zeros(K-1,1);
m2 = zeros(K-1,1);
while (k2 < K)
    m1(k2) = mk(1);
    m2(k2) = mk(2); 
    beta = (k2-1)/(k2+2);
    mk = nk_1 - s3 * (A * nk_1 + (b + sign(mk)));
    nk = mk + beta * (mk - mk_1);        
    mk_1 = mk;
    nk_1 = nk;         
    k2 = k2 + 1;
end

figure(1)
plot(x1,x2,'-.','color','red')
hold on
plot(x1(ceil(1/sqrt(s1))),x2(ceil(1/sqrt(s1))),'ro')
hold on
plot(x1(ceil(2/sqrt(s1))),x2(ceil(2/sqrt(s1))),'r+')
hold on
plot(x1(ceil(3/sqrt(s1))),x2(ceil(3/sqrt(s1))),'r^')
hold on
plot(u1,u2,'-.','color','blue')
hold on
plot(u1(ceil(1/sqrt(s2))),u2(ceil(1/sqrt(s2))),'bo')
hold on
plot(u1(ceil(2/sqrt(s2))),u2(ceil(2/sqrt(s2))),'b+')
hold on
plot(u1(ceil(3/sqrt(s2))),u2(ceil(3/sqrt(s2))),'b^')
hold on
plot(m1,m2,'-.','color','black')
hold on
plot(m1(ceil(1/sqrt(s3))),m2(ceil(1/sqrt(s3))),'o','color','black')
hold on
plot(m1(ceil(2/sqrt(s3))),m2(ceil(2/sqrt(s3))),'+','color','black')
hold on
plot(m1(ceil(3/sqrt(s3))),m2(ceil(3/sqrt(s3))),'^','color','black')
xlabel('x1')
ylabel('x2')

test = a' * ( a*yk - [4 2 0]' ) + [-1 1]';