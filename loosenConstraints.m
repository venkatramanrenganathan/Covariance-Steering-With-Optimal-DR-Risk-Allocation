function PS = loosenConstraints(PS,Nactive,deltaRes,bool)
    N = PS.N;
    Ns = PS.Ns;
    deltax = PS.deltax;
    
    for k = 1:N
        for j = 1:Ns
            if bool(j,k)
                deltax(j,k) = deltax(j,k) + deltaRes / Nactive;
            end
        end
    end
    
    PS.deltax = deltax;
end