function deltaTrue = computeTrueDeltaCC(PS,V,K,k)
    nx = PS.nx;
    N = PS.N;
    ScriptA = PS.ScriptA;
    ScriptB = PS.ScriptB;
    ScriptD = PS.ScriptD;
    mu0 = PS.mu0;
    Sigma0 = PS.Sigma0;
    offset = PS.stateCC_offset;
    
    SigmaY = ScriptA*Sigma0*ScriptA'+ScriptD*ScriptD';
    IplusBK = eye((N+1)*nx) + ScriptB * K;
    SigmaX = IplusBK * SigmaY * IplusBK';
    EX = ScriptA * mu0 + ScriptB * V;
    Ek = [zeros(nx,k*nx) eye(nx) zeros(nx,(N-k)*nx)];
    SigmaXk = Ek * SigmaX * Ek';
    
    a = PS.a;
    c = PS.c;
    beta = PS.beta;
    gamma = PS.gamma;
    
    EXk = Ek * EX;
    rk = gamma * EXk(2) + offset;
    fk = rk / sqrt(2);

    inside1 = (-fk - a' * EXk) / sqrt(a' * SigmaXk * a);
    inside2 = (fk - a' * EXk) / sqrt(a' * SigmaXk * a);
    inside3 = (-fk - c' * EXk) / sqrt(a' * SigmaXk * a);
    inside4 = (fk - c' * EXk) / sqrt(a' * SigmaXk * a);
    
    deltaTrue1 = (normcdf(inside1)) / beta;
    deltaTrue2 = (1 - normcdf(inside2)) / beta;
    deltaTrue3 = (normcdf(inside3)) / (1 - beta);
    deltaTrue4 = (1 - normcdf(inside4)) / (1 - beta);
    
    deltaTrue = [deltaTrue1;deltaTrue2;deltaTrue3;deltaTrue4];
end