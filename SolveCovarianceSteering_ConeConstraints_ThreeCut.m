%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file solves the OCS problem with state chance constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Kazuhide Okamoto, DCSL, 2019
% (c) Joshua Pilipovsky, DCSL, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [time,obj,Vsol,Ksol,fsol,tsol] = SolveCovarianceSteering_ConeConstraints_ThreeCut(PS)   
    % % % Lower Stage Optimization % % %
    
    % System parameters
    N = PS.N;
    nx = PS.nx;
    nu = PS.nu;
    ScriptA = PS.ScriptA;
    ScriptB = PS.ScriptB;
    ScriptD = PS.ScriptD;
    mu0 = PS.mu0;
    muf = PS.muf;
    Sigma0 = PS.Sigma0;
    Sigmaf = PS.Sigmaf;
    Q = PS.Q;
    R = PS.R;
    
    % Define optimization variables
    V = sdpvar(nu*N,1);
    K = [];
    for i = 1:N
        K = blkdiag(K,sdpvar(nu,nx));
    end
    K = [K,zeros(nu*N,nx)];
    f = sdpvar(nx,N);
    t = sdpvar(nx,N);
    
    % Useful matrices
    SigmaY = ScriptA * Sigma0 * ScriptA' + ScriptD * ScriptD';
    IplusBK = eye((N+1)*nx) + ScriptB*K;
    EN = [zeros(nx,N*nx) eye(nx)];
    S = computeS(SigmaY);

    % System Dynamics
    EX = ScriptA * mu0 + ScriptB * V;
    Z = EN * IplusBK * S;

    % Define constraints
    Constraints = [];
    
    % Boundary condition
    Constraints = [Constraints, EN * EX == muf];      
    Constraints = [Constraints, [Sigmaf Z; Z' eye((N+1)*nx)] >= 0];
    
    % Chance Constraints
    Constraints = chanceConstraints(Constraints,EX,S,IplusBK,V,K,f,t,PS); 
    
    % Objective Function
    Jmean = EX' * Q * EX + V' * R * V;
    Jcov = trace((IplusBK' * Q * IplusBK + K' * R * K) * SigmaY);
    Objective = Jmean + Jcov;
    
    % Solve the Problem
    options = sdpsettings('solver','mosek');
    sol = optimize(Constraints,Objective,options);
    
    % Get performance measure values
    time = sol.solvertime;
    obj = value(Objective);
    
    % Get values of optimization variables
    Vsol = value(V);
    Ksol = value(K);
    fsol = value(f);
    tsol = value(t);
end

function Constraints = chanceConstraints(Constraints,EX,S,IplusBK,V,K,f,t,PS)
    % Control constraint parameters
    a_u = PS.alphau;
    b_u = PS.betau;
    Nc = PS.Nc;
    N = PS.N;
    deltau = PS.deltau;
    
    % System dimensions
    nx = PS.nx;
    nu = PS.nu;
    
    % State constraint parameters
    beta = PS.beta;
    A = PS.A_cone;
    b = PS.b;
    c = PS.c;
    d = PS.d;    
    deltax = PS.deltax;
    
    
    % % % Control chance constraints % % %
    
    for k = 0:N-1
        Fk = [zeros(nu,k*nu) eye(nu) zeros(nu,(N-k-1)*nu)];
        for i = 1:Nc
            ICDFu = norminv(1 - deltau(i,k+1));
            controlConstraint = a_u(:,i)' * Fk * V + ICDFu * norm(S' * K' * Fk' * a_u(:,i)) <= b_u(i);
            Constraints = [Constraints,controlConstraint];
        end
    end
    
    % % % State chance constraints % % %
    
    for k = 1:N
        Ek = [zeros(nx,k*nx) eye(nx) zeros(nx,(N-k)*nx)];
        for i = 1:nx
            stateConstraint1 = t(i,k) >= norm(S' * IplusBK' * Ek' * A(i,:)');
            stateConstraint2 = -(f(i,k) + b(i) + A(i,:) * Ek * EX) <= norminv(beta * deltax(k)) * t(i,k);
            stateConstraint3 = f(i,k) - b(i) - A(i,:) * Ek * EX >= norminv(1 - beta * deltax(k)) * t(i,k);
            stateConstraint4 = -f(i,k) <= norminv(beta * deltax(k) / 2) * t(i,k);
            stateConstraints = [stateConstraint1,stateConstraint2,stateConstraint3,stateConstraint4];
            Constraints = [Constraints,stateConstraints];
        end
        fk = f(:,k);
        Constraints = [Constraints, norm(fk) <= c' * Ek * EX + d];
    end

end