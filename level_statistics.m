% gamma is irrational, and approximated by F(m-1)/F(m)
% system length is L = F(m)
% fibonacci sequence
m = 15;
L = fibonacci(m);
gamma = fibonacci(m-1)/fibonacci(m);
n = 1:L;

% onsite quasiperiodic potential
V = 400;
t1 = ones(L-1,1);
en = zeros(L,300);
for ind = 1:300
    %V1 = (rand(L,1)-0.5)*50;
    phi = rand;
    V1 = V*cos(2*pi*(gamma*n+phi));
    H = diag(V1) + diag(t1,1) + diag(t1,-1);
    % periodic boundary condition
    H(1,L) = 1;
    H(L,1) = 1;
    % get the true ground energy and true ground state
    [v,d] = eig(H);
    Eval = diag(d);
    en(:,ind) = Eval;
end
E = en(:);
histogram(E,500)