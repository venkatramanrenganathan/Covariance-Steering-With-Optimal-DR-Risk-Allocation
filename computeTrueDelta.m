function deltaTrue = computeTrueDelta(PS,V,K,j,k)
    a_x = PS.alphax;
    b_x = PS.betax;
    nx = PS.nx;
    N = PS.N;
    ScriptA = PS.ScriptA;
    ScriptB = PS.ScriptB;
    ScriptD = PS.ScriptD;
    mu0 = PS.mu0;
    Sigma0 = PS.Sigma0;
    
    SigmaY = ScriptA*Sigma0*ScriptA'+ScriptD*ScriptD';
    IplusBK = eye((N+1)*nx) + ScriptB * K;
    SigmaX = IplusBK * SigmaY * IplusBK';
    EX = ScriptA * mu0 + ScriptB * V;
    Ek = [zeros(nx,k*nx) eye(nx) zeros(nx,(N-k)*nx)];

    inside = (b_x(j) - a_x(:,j)' * Ek * EX) / sqrt(a_x(:,j)' * Ek * SigmaX * Ek' * a_x(:,j));
    deltaTrue = 1 - normcdf(inside);
end