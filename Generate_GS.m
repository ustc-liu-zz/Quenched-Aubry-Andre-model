%% propagate i*C'= H*C in imag-time to get the ground state 
% prepare the Hamiltonian matrix
t = 1;
% gamma is irrational, and approximated by F(m-1)/F(m)
% system length is L = F(m)
% fibonacci sequence
m = 16;
L = fibonacci(m);
gamma = fibonacci(m-1)/fibonacci(m);
phi = rand(1);
n = 1:L;

% onsite quasiperiodic potential
V = 2.5*t;
V1 = V*cos(2*pi*(gamma*n+phi));
t1 = t*ones(L-1,1);
H = diag(V1) + diag(t1,1) + diag(t1,-1);
% periodic boundary condition
H(1,L) = t;
H(L,1) = t;
% get the true ground energy and true ground state
[v,d] = eig(H);

% prepare the initial wave-function
%C = zeros(L,1);
%C(floor(L/2)) = 1;
C = ones(L,1)/sqrt(L);

% time step
h = 0.01;
Eval_b = 0;
Eval_a = C'*H*C;
tol = 1e-12;
% fourth-order Runge-Kutta
while abs(Eval_b - Eval_a) >= tol
    Eval_b = Eval_a;
    k1 = -h*H*C;
    k2 = -h*H*(C+k1/2);
    k3 = -h*H*(C+k2/2);
    k4 = -h*H*(C+k3);
    C = C + k1/6 + k2/3 + k3/3 + k4/6;
    NN = sum(abs(C).^2);
    C = C/sqrt(NN);
    Eval_a = C'*H*C;
end
plot(n,abs(C).^2,n,abs(v(:,1)).^2,'-.','linewidth',1.5)
save('GS.mat','t','m','phi','C')