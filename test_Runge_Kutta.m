%% test fourth-order Runge-Kutta method
% solve y' = y, error is O(h^5)
% initial condition
N = 10;
y = zeros(N,1);
y(1) = 1;
h = 0.1;
% fourth-order Runge-Kutta
for ind = 2:N
    k1 = h*y(ind-1);
    k2 = h*(y(ind-1)+k1/2);
    k3 = h*(y(ind-1)+k2/2);
    k4 = h*(y(ind-1)+k3);
    y(ind) = y(ind-1) + k1/6 + k2/3 + k3/3 + k4/6;
end
x = 0:h:(N*h-h);
% exact solution
z = exp(x);
plot(x,y,'o',x,z)