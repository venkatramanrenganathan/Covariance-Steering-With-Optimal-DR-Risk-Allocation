function deltaTrue = computeTrueDeltaC_ConeConstraints_RUB(PS,V,K,f)    
% System parameters
ScriptA = PS.ScriptA;
ScriptB = PS.ScriptB;
ScriptD = PS.ScriptD;

mu0 = PS.mu0;
Sigma0 = PS.Sigma0;

N = PS.N;
nx = PS.nx;
b = PS.b;
beta = PS.beta;
A = PS.A_cone;

% Useful matrices
SigmaY = ScriptA * Sigma0 * ScriptA' + ScriptD * ScriptD';
IplusBK = eye((N+1)*nx) + ScriptB*K;
S = computeS(SigmaY);

% Nominal state trajectory
EX = ScriptA * mu0 + ScriptB * V;

% Compute true risks
deltaTrue = zeros(2,nx,N);

for k = 1:N
    Ek = [zeros(nx,k*nx) eye(nx) zeros(nx,(N-k)*nx)];
    for i = 1:nx
        normterm1 = (f(i,k) - b(i) - A(i,:) * Ek * EX) / norm(S' * IplusBK' * Ek' * A(i,:)');
        normterm2 = (f(i,k) + b(i) + A(i,:) * Ek * EX) / norm(S' * IplusBK' * Ek' * A(i,:)');
        deltaTrue1 = (2 / beta) * (1 - normcdf(normterm1));
        deltaTrue2 = (2 / beta) * (1 - normcdf(normterm2));
        deltaTruei = [deltaTrue1; deltaTrue2];
        deltaTrue(:,i,k) = deltaTruei;
    end
end

end