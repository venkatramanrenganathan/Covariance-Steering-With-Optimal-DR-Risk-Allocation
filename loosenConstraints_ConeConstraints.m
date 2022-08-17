function PS = loosenConstraints_ConeConstraints(PS,Nactive,deltaRes,bool)
N = PS.N;
deltax = PS.deltax;

for k = 1:N
    if bool(k)
        deltax(k) = deltax(k) + deltaRes / Nactive;
    end
end

PS.deltax = deltax;
end