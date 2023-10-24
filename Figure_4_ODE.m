clc
clear
x0 = [1,0]';
tspan = [0.01 9];
options = odeset('RelTol',1e-12,'Stats','on','OutputFcn',@odeplot);
[t,x] = ode45(@(t,x) odefcn1(t,x), tspan, x0, options);

error_ = zeros(length(x),1);
x_ = x( :, 1);

k = 1;
while (k < length(t))
    error_(k) = t(k)^2 * abs(x_(k)); 
    k = k+1;
end

figure(1)
plot(t',error_)
xlabel('t')
ylabel('t^2(f - f*)')

function dxdt = odefcn1(t,x)
  r = 4;
  dxdt = ones(2,1);
  dxdt(1) = x(2);
  
  if  (x(1) < 0)
    dxdt(2) = (-r/t) * dxdt(1) + 1;
    
  elseif (x(1) > 0)
    dxdt(2) = (-r/t) * dxdt(1) - 1;
    
  else
    dxdt(2) = (-r/t) * dxdt(1) ;
  
  end  
  
end
