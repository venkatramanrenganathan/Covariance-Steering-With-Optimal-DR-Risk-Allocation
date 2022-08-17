function Debugging(PS)
ScriptA = PS.ScriptA;
ScriptB = PS.ScriptB;
ScriptD = PS.ScriptD;
V = PS.V;
K = PS.K;
mu0 = PS.mu0;
Sigma0 = PS.Sigma0;
N = PS.N;
nx = PS.nx;
t = PS.t;
f = PS.f;
SigmaY = ScriptA * Sigma0 * ScriptA' + ScriptD * ScriptD';
S = computeS(SigmaY);
b = PS.b;
beta = PS.beta;
deltax = PS.deltax;
A = PS.A_cone;
c = PS.c;
d = PS.d;

EX = ScriptA * mu0 + ScriptB * V;
IplusBK = eye((N+1)*nx) + ScriptB*K;
Constraints = [];
for k = 1:N
    Ek = [zeros(nx,k*nx) eye(nx) zeros(nx,(N-k)*nx)];
    for i = 1:nx
        stateConstraint1 = t - norm(S' * IplusBK' * Ek' * A(i,:)'); 
        stateConstraint2 = norminv(beta * deltax(k)) * t + f(i) + b(i) + A(i,:) * Ek * EX; 
        stateConstraint3 = f(i) - b(i) - A(i,:) * Ek * EX - norminv(1 - beta * deltax(k)) * t; 
        stateConstraint4 = f(i) + norminv(beta * deltax(k) / 2) * t;
        stateConstraints = [stateConstraint1,stateConstraint2,stateConstraint3,stateConstraint4];
        Constraints = [Constraints,stateConstraints];
    end
    Constraints = [Constraints, c' * Ek * EX + d - norm(f)];
end

end