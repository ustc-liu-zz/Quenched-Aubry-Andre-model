S = zeros(Nt/100,1);
n = n';
for ind = 1:Nt/100
    S(ind) = sqrt(sum(n.^2.*abs(Evec(:,ind)).^2)-(sum(n.*abs(Evec(:,ind)).^2))^2)/L;
end