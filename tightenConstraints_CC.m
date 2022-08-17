function PS = tightenConstraints_CC(PS,alpha,V,K,bool,char)
    N = PS.N;
    deltax = PS.deltax;
    
    if strcmp(char,'GEO')
        for k = 1:N
            if ~bool(k)
                deltaTruek = computeTrueDeltaCC2(PS,V,K,k);
                deltax(k) = alpha * deltax(k) + (1 - alpha) * deltaTruek;
            end
        end
    else
        for k = 1:N
            if ~bool(k)
                deltaTruek = max(computeTrueDeltaCC(PS,V,K,k));
                deltax(k) = alpha * deltax(k) + (1 - alpha) * deltaTruek;
            end
        end
    end
    
    PS.deltax = deltax;
end