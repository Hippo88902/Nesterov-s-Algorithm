x0 = 1;
K = 10000;
s = 10^-6;
x_k1 = x0;
y_k1 = x0;
xk = zeros(K-1,1);
t = ones(K-1,1);
error = zeros(K-1,1);

r = 4; %r值
k = 1;
while (k < K)
    
    beta = (k-1)/(k+r-1);
    if(y_k1 > 0)
        x_k = y_k1 - s ;
    elseif(y_k1 < 0)
        x_k = y_k1 + s ;
    else
        x_k = y_k1;
    end    
    y_k = x_k + beta*(x_k - x_k1);
    xk(k) = x_k;
    t(k) = k*sqrt(s);
    error(k) = t(k)^2 * abs(x_k);
    y_k1 = y_k;
    x_k1 = x_k;     
    k = k + 1;

end

ki = linspace(1,K-1,K-1);
figure(1)
plot(ki,error','-.','color','red')
xlabel('iteration')
ylabel('sk^2(f - f*)')


