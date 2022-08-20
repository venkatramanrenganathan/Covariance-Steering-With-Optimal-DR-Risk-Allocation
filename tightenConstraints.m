function PS = tightenConstraints(PS,alpha,V,K,bool, riskSelectFlag)
    N = PS.N;
    Ns = PS.Ns;
    deltax = PS.deltax;
    
    for k = 1:N
        for j = 1:Ns
            if ~bool(j,k)
                if(riskSelectFlag == 1)
                    % Gaussian Case
                    deltaTrue = computeTrueDelta(PS,V,K,j,k);
                elseif(riskSelectFlag == 2)
                    % DR Case
                    deltaTrue = computeDRTrueDelta(PS,V,K,j,k);
                end                                
                deltax(j,k) = alpha * deltax(j,k) + (1 - alpha) * deltaTrue;
            end
        end
    end
    
    PS.deltax = deltax;
end