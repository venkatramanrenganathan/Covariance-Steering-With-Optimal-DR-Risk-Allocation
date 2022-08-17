function deltaTrue = computeTrueDeltaC_ConeConstraints_ThreeCut(PS,V,f,t)    
ScriptA = PS.ScriptA;
ScriptB = PS.ScriptB;
mu0 = PS.mu0;
N = PS.N;
nx = PS.nx;
b = PS.b;
beta = PS.beta;
A = PS.A_cone;

EX = ScriptA * mu0 + ScriptB * V;

deltaTrue = zeros(3,nx,N);

for k = 1:N
    Ek = [zeros(nx,k*nx) eye(nx) zeros(nx,(N-k)*nx)];
    for i = 1:nx
        normterm1 = -(1 / t(i,k)) * (f(i,k) + b(i) + A(i,:) * Ek * EX);
        normterm2 = (1 / t(i,k)) * (f(i,k) - b(i) - A(i,:) * Ek * EX);
        normterm3 = -f(i,k) / t(i,k);
        deltaTrue1 = (1 / beta) * normcdf(normterm1);
        deltaTrue2 = (1 / beta) * (1 - normcdf(normterm2));
        deltaTrue3 = (2 / beta) * normcdf(normterm3);
        deltaTruei = [deltaTrue1; deltaTrue2; deltaTrue3];
        deltaTrue(:,i,k) = deltaTruei;
    end
end

end