%% test fourth-order Predictor-Corrector Runge-Kutta method
% solve y' = y, error is O(h^5)
% initial condition
N = 200;
y = zeros(N,1);
y(1) = 1;
h = 0.01;
% fourth-order Runge-Kutta
for ind = 2:4
    k1 = h*y(ind-1);
    k2 = h*(y(ind-1)+k1/2);
    k3 = h*(y(ind-1)+k2/2);
    k4 = h*(y(ind-1)+k3);
    y(ind) = y(ind-1) + k1/6 + k2/3 + k3/3 + k4/6;
end

for ind = 5:N
    temp1 = 55*y(ind-1)-59*y(ind-2)+37*y(ind-3);
    temp2 = -9*y(ind-4);
    W0 = y(ind-1) + h*(temp1+temp2)/24;
    temp3 = 9*W0+19*y(ind-1)-5*y(ind-2)+y(ind-3);
    y(ind) = y(ind-1) + h*temp3/24;
end
x = 0:h:(N*h-h);
% exact solution
z = exp(x);
plot(x,y,'o',x,z)