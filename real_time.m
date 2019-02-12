%% real time quench process
m = 13;
L = fibonacci(m);
gamma = fibonacci(m-1)/fibonacci(m);
n = 1:L;
n = n';

t = 1;
V = 3;
phi = rand;
V1 = V*cos(2*pi*(gamma*n+phi));
t1 = t*ones(L-1,1);
H = diag(V1) + diag(t1,1) + diag(t1,-1);
% periodic boundary condition
H(1,L) = t;
H(L,1) = t;
[v,d] = eig(H);
% initial state
C = v(:,1);
Cg = C;
% quench time
tauq = 64;

% imaginary evolution time step
dt = 0.01*1i;

% evolution steps
Nt = 10*tauq/(dt*(-1i));

% steps to record data
Nd = 10;

% store energy
Eval = zeros(Nt,1);

% wave pocket width
Delta_w = zeros(Nt/Nd,1);

% store wave-function
Evec = zeros(L,Nt/Nd);
indy = 1;
for ind = 1:Nt
    V = -tanh((ind*dt*(-1i)-5*tauq)/tauq)+2;
    V1 = V*cos(2*pi*(gamma*n+phi));
    t1 = t*ones(L-1,1);
    H = diag(V1) + diag(t1,1) + diag(t1,-1);
    % periodic boundary condition
    H(1,L) = t;
    H(L,1) = t;
    
    % fourth-order Runge-Kutta
    k1 = -dt*H*C;
    k2 = -dt*H*(C+k1/2);
    k3 = -dt*H*(C+k2/2);
    k4 = -dt*H*(C+k3);
    C = C + k1/6 + k2/3 + k3/3 + k4/6;
    Eval(ind) = C'*H*C;
    if mod(ind,Nd) == 0
        Evec(:,indy) = C;
        Delta_w(indy) = sqrt(sum(n.^2.*abs(C).^2)-sum(n.*abs(C).^2)^2);
        indy = indy + 1;
    end
end