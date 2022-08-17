function PS = loosenConstraints_CC(PS,Nactive,deltaRes,bool,char)
    N = PS.N;
    deltax = PS.deltax;
    
    if strcmp(char,'GEO')
        for k = 1:N
            if bool(k)
                deltax(k) = deltax(k) + deltaRes / Nactive;
            end
        end
    else
        for k = 1:N
            if bool(k)
                deltax(k) = deltax(k) + deltaRes / Nactive;
            end
        end
    end
    
    PS.deltax = deltax;
end