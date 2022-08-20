function [bool,diff] = isActive(V,K,PS, riskSelectFlag)
    eps = 1E-06;
    Ns = PS.Ns;
    N = PS.N;
    bool = false(Ns,N);
    
    a_x = PS.alphax;
    b_x = PS.betax;
    deltax = PS.deltax;
    nx = PS.nx;
    ScriptA = PS.ScriptA;
    ScriptB = PS.ScriptB;
    ScriptD = PS.ScriptD;
    
    mu0 = PS.mu0;
    Sigma0 = PS.Sigma0;
    SigmaY = ScriptA*Sigma0*ScriptA'+ScriptD*ScriptD';
    S = computeS(SigmaY);
    
    EX = ScriptA * mu0 + ScriptB * V;
    IplusBK = eye((N+1)*nx) + ScriptB * K;
    
    for k = 1:N
        Ek = [zeros(nx,k*nx) eye(nx) zeros(nx,(N-k)*nx)];
        for j = 1:Ns
            if(riskSelectFlag == 1)
                % Gaussian Case
                ICDF = norminv(1-deltax(j,k));
            elseif(riskSelectFlag == 2)
                % DR Case
                ICDF = sqrt((1-deltax(j,k))/(deltax(j,k)));
            end
            diff(j,k) = b_x(j) - (a_x(:,j)'*Ek*EX + ICDF*norm(S'*IplusBK'*Ek'*a_x(:,j)));
            if  diff(j,k) < eps
                bool(j,k) = true;
            end
        end
    end
end