function [bool,diff] = isActive_CC(V,K,PS,char)
    eps = 1E-06;
    N = PS.N;
    if strcmp(char,'GEO')
        bool = false(1,N);
    else
        bool = false(4,N);
    end

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
    
    Ac = PS.A_cone(1:3,1:3);
    cc = PS.c(1:3);
    H = [1 0 0; 0 0 1];
    Ip = [eye(3), zeros(3)];
    d = PS.d;
    
    if strcmp(char,'GEO')
        for k = 1:N
            Ek = [zeros(nx,k*nx) eye(nx) zeros(nx,(N-k)*nx)];
            EXk = Ek * EX;        
            xi_bar_k = H * Ac * Ip * Ek * EX;
            Rk = (cc' * Ip * EXk + d) - norm(xi_bar_k);
            diff(k) = Rk - sqrt(2 * log(1 / deltax(k))) * norm(S' * IplusBK' * Ek' * Ip' * Ac' * H');

            if  diff(k) < eps
                bool(k) = true;
            end
        end
    else
        for k = 1:N
            Ek = [zeros(nx,k*nx) eye(nx) zeros(nx,(N-k)*nx)];
            EXk = Ek * EX;
            rk = gamma * EXk(2) + offset;
            fk = rk / sqrt(2);
            inside1 = (-fk - a' * EXk) / norm(S' * IplusBK' * Ek' * a);
            inside2 = (fk - a' * EXk) / norm(S' * IplusBK' * Ek' * a);
            inside3 = (-fk - c' * EXk) / norm(S' * IplusBK' * Ek' * a);
            inside4 = (fk - c' * EXk) / norm(S' * IplusBK' * Ek' * a);

            term1 = (1 / beta) * normcdf(inside1);
            term2 = (1 / beta) * (1 - normcdf(inside2));
            term3 = (1 / (1 - beta)) * normcdf(inside3);
            term4 = (1 / (1 - beta)) * (1 - normcdf(inside4));

            diff(1,k) = deltax(k) - term1;
            diff(2,k) = deltax(k) - term2;
            diff(3,k) = deltax(k) - term3;
            diff(4,k) = deltax(k) - term4;
            
            for j = 1:4
                if  diff(j,k) < eps
                    bool(j,k) = true;
                end
            end
        end
        
        bool = any(bool);
    end
end