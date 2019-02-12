%% real time quench process
% real evolution time step
dt = 0.01*1i;
t = 1;
m = 15;
L = fibonacci(m);
gamma = fibonacci(m-1)/fibonacci(m);
n = 1:L;
n = n';

% initial V
V = 3;
t1 = t*ones(L-1,1);
% steps to record data
Nd = 10;
% repeat times
Nr = 5;
xi = zeros(6,Nr);
indz = 1;
for tauq = [16 32 64 128 256 512]
    for indr = 1:Nr        
        phi = rand;
        V1 = V*cos(2*pi*(gamma*n+phi));        
        H = diag(V1) + diag(t1,1) + diag(t1,-1);
        % periodic boundary condition
        H(1,L) = t;
        H(L,1) = t;
        [v,d] = eig(H);
        % initial state
        C = v(:,1);
        
        % evolution steps
        Nt = 10*tauq/(dt*(-1i));
        
        % wave pocket width
        Delta_w = zeros(Nt/Nd,1);
        indy = 1;
        for ind = 0:Nt-1
            % quench parameter V(t) = -tanh(t/(tauq))+2, t0 = -5*tauq
            V_temp1 = -tanh((ind*dt*(-1i)-5*tauq)/tauq)+2;
            V_temp2 = -tanh(((ind+0.5)*dt*(-1i)-5*tauq)/tauq)+2;
            V_temp3 = -tanh(((ind+1)*dt*(-1i)-5*tauq)/tauq)+2;
            V1 = V_temp1*cos(2*pi*(gamma*n+phi));
            V2 = V_temp2*cos(2*pi*(gamma*n+phi));
            V3 = V_temp3*cos(2*pi*(gamma*n+phi));
            t1 = t*ones(L-1,1);
            H1 = diag(V1) + diag(t1,1) + diag(t1,-1);
            H2 = diag(V2) + diag(t1,1) + diag(t1,-1);
            H3 = diag(V3) + diag(t1,1) + diag(t1,-1);
            % periodic boundary condition
            H1(1,L) = t;
            H1(L,1) = t;
            
            H2(1,L) = t;
            H2(L,1) = t;
            
            H3(1,L) = t;
            H3(L,1) = t;
            
            % fourth-order Runge-Kutta
            k1 = -dt*H1*C;
            k2 = -dt*H2*(C+k1/2);
            k3 = -dt*H2*(C+k2/2);
            k4 = -dt*H3*(C+k3);
            C = C + k1/6 + k2/3 + k3/3 + k4/6;            
            if mod(ind,Nd) == 0
                Delta_w(indy) = sqrt(sum(n.^2.*abs(C).^2)-sum(n.*abs(C).^2)^2);
                indy = indy + 1;
            end
        end
        xi(indz,indr) = Delta_w(Nt/(2*Nd));
    end
    indz = indz + 1
end
tauq = [16 32 64 128 256 512]';
save('xi_tauq.mat','tauq','xi')