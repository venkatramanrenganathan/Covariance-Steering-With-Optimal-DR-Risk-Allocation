function deltaTrue = computeTrueDeltaCC2(PS,V,K,k)

    nx = PS.nx;
    N = PS.N;
    
    ScriptA = PS.ScriptA;
    ScriptB = PS.ScriptB;
    ScriptD = PS.ScriptD;
    
    mu0 = PS.mu0;
    Sigma0 = PS.Sigma0;
    
    SigmaY = ScriptA*Sigma0*ScriptA'+ScriptD*ScriptD';
    S = computeS(SigmaY);
    IplusBK = eye((N+1)*nx) + ScriptB * K;
    
    Ek = [zeros(nx,k*nx) eye(nx) zeros(nx,(N-k)*nx)];
    EX = ScriptA * mu0 + ScriptB * V;
    EXk = Ek * EX;

    Ac = PS.A_cone(1:3,1:3);
    cc = PS.c(1:3);
    H = [1 0 0; 0 0 1];
    Ip = [eye(3), zeros(3)];
    d = PS.d;
    
    xi_bar_k = H * Ac * Ip * Ek * EX;
    Rk = (cc' * Ip * EXk + d) - norm(xi_bar_k);
    
    normSquared = norm(S' * IplusBK' * Ek' * Ip' * Ac' * H')^2;
    
    deltaTrue = exp(-0.5 * Rk^2 / normSquared);
end