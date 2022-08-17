function volume = computeVolumeN(PS,K)
    alpha = 1E05;
    N = PS.N;
    nx = PS.nx;
    SigmaY = PS.ScriptA*PS.Sigma0*PS.ScriptA'+PS.ScriptD*PS.ScriptD';
    IplusBK = eye((N+1)*nx)+PS.ScriptB*K;
    SigmaX = IplusBK * SigmaY * IplusBK';
    
    EN = [zeros(nx,N*nx) eye(nx)];
    SigmaN = EN*SigmaX*EN';
    
    volume = det(SigmaN(1:3,1:3));
    %volume = alpha * det(SigmaN);
end