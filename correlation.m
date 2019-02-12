Np = 987;
nr = 1;   % number of random realizations of phi
t = 40;

xii = zeros(1,t);

for k = 1:t
    
    W = 1 + (0.3/(1.09).^(k-1)); % epsilon = (0.3/(1.03).^k); k runs from 1 to 50
    
    
    for fg = 1:nr
        
        %phi = 0;
        phi = random('unif',0,1);
        H=AA_Ham(Np,W,phi); %see the local function generating the tight binding matrix
        %of AA model at desired potential strength
        [Wf,D]=eig(H); %find eigenvalues and eigenstates of tb matrix
        
        gss = Wf(:,1); % it gives the groundstate at that disorder strength
        
        gsf = abs(gss).^2;
        
        means = 0;
        
        for kk = 1:Np
            mmss = kk*gsf(kk);
            means = means + mmss;
        end
        
        gsf = circshift(gsf,round(Np/2 - means));
        
        mean=0;
        
        for lk = 1:Np
            mms = lk*gsf(lk);
            mean = mean + mms;
        end
        
        var = 0;
        
        for l = 1:Np
            knk = ((l-mean).^2)*gsf(l);
            var = var + knk;
        end
        std = sqrt(var);
        
        xii(k) = xii(k) + std;
    end
    
end

xi = xii/nr;

esp = zeros(1,t);
for l = 1:t
    esp(l) = (0.3/(1.09).^(l-1));
end

plot(log10(esp),log10(xi),'g*')

fileID = fopen('correlation.txt','w');

for ty = 1:t
    fprintf(fileID,'%f %f\n',esp(ty),xi(ty));
end
fclose(fileID);






function y = AA_Ham(Np,W,phi) %Setting up the tb matrix from which we calculate
%the initial localized groundstate
t=1; %Nearest-neighbor hopping as the unit of energy
irr= 610/987;
%irr = (sqrt(5)-1)/2;

H =-(t*diag(ones(1,Np-1),1))-(t*diag(ones(1,Np-1),-1));
H(1,Np)=-t; H(Np,1)=-t;

U=2*W*cos(2*pi*((1:Np)*irr + phi));


H = H +diag(U);
y=H;
end

