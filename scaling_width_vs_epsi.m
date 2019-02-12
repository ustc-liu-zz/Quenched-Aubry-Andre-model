%% scaling of wave pocket width vs distance \espilison=V-2t
t = 1;
m = 16;
L = fibonacci(m);
L_cut = floor(L/2);
gamma = fibonacci(m-1)/fibonacci(m);
n = 1:L;
n = n';

x = logspace(-1.5,-0.3,30)+2;
NN = length(x);
indy = 1;
ft = fittype('a*x+b');
phi = rand;
t1 = t*ones(L-1,1);
% wave pocket width
Delta_w = zeros(NN,1);
for V = x
    V1 = V*cos(2*pi*(gamma*n+phi));
    H = diag(V1) + diag(t1,1) + diag(t1,-1);
    H(1,L) = t;
    H(L,1) = t;
    [Evec,d] = eig(H);
    C = Evec(:,1);
    Delta_w(indy) = sqrt(sum(n.^2.*abs(C).^2)-sum(n.*abs(C).^2)^2);
    indy = indy + 1;
    
end
epsilon = x-2;
epsilon = epsilon';

y = 0.01:0.01:1;
myfit = fit(log(epsilon),log(Delta_w),ft,'StartPoint',[-1 1]);
z = y.^(myfit.a)*exp(myfit.b);
loglog(epsilon,Delta_w,'d',y,z,'k-')
mytitle = join(['wave pocket width: \xi\sim \epsilon^{-\nu}, \nu =  ',num2str(myfit.a)]);
xlabel('\epsilon')
ylabel('\xi')
title(mytitle)
xlim([0.01 0.8])
ylim([0.5 50])