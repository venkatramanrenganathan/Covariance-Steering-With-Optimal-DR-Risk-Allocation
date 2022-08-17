function PS = tightenConstraints(PS,alpha,V,K,bool)
    N = PS.N;
    Ns = PS.Ns;
    deltax = PS.deltax;
    
    for k = 1:N
        for j = 1:Ns
            if ~bool(j,k)
                deltaTrue = computeTrueDelta(PS,V,K,j,k);
                deltax(j,k) = alpha * deltax(j,k) + (1 - alpha) * deltaTrue;
            end
        end
    end
    
    PS.deltax = deltax;
end