%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file solves the OCS problem with state chance constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Kazuhide Okamoto, DCSL, 2019
% (C) Joshua Pilipovsky, DCSL, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [time,obj,Vsol,Ksol] = SolveCovarianceSteeringCC(PS,char)   
    %%% Lower Stage %%%
    N = PS.N;
    nx = PS.nx;
    nu = PS.nu;
    ScriptA = PS.ScriptA;
    ScriptB = PS.ScriptB;
    ScriptD = PS.ScriptD;
    
    % Define optimization variables
    V = sdpvar(nu*N,1);
    K = [];
    for i = 1:N
        K = blkdiag(K,sdpvar(nu,nx));
    end
    K = [K,zeros(nu*N,nx)];
    
    % Useful matrices
    SigmaY = ScriptA*PS.Sigma0*ScriptA'+ScriptD*ScriptD';
    IplusBK = eye((N+1)*nx)+ScriptB*K;
    EN = [zeros(nx,N*nx) eye(nx)];
    S = computeS(SigmaY);

    % System Dynamics
    EX = ScriptA * PS.mu0 + ScriptB * V;
    Z = EN*IplusBK*S;

    % Define constraints
    Constraints = [];
    
    % Boundary condition
    Constraints = [Constraints, EN*EX == PS.muf];      
    Constraints = [Constraints, [PS.Sigmaf Z; Z' eye((N+1)*nx)] >= 0];
    
    % Chance Constraints
    Constraints = chanceConstraints(Constraints,EX,S,SigmaY,IplusBK,V,K,PS,char); 
    
    % Objective Function
    Jmean = EX'*PS.Q*EX + V'*PS.R*V;
    Jcov = trace((IplusBK'*PS.Q*IplusBK + K'*PS.R*K)*SigmaY);
    Objective = Jmean + Jcov;
    
    % Solve the Problem
    options = sdpsettings('solver','mosek');
    sol = optimize(Constraints,Objective,options);
    
    % Get performance measure values
    time = sol.solvertime;
    obj = value(Objective);
    
    % Get control values
    Vsol = value(V);
    Ksol = value(K);
end

function Constraints = chanceConstraints(Constraints,EX,S,SigmaY,IplusBK,V,K,PS,char)
    a_u = PS.alphau;
    b_u = PS.betau;
    Nc = PS.Nc;
    N = PS.N;
    deltax = PS.deltax;
    deltau = PS.deltau;
    nx = PS.nx;
    nu = PS.nu;
    
    beta = PS.beta;
    A = PS.A_cone;
    b = PS.b;
    c = PS.c;
    d = PS.d;    
    
    % % % Control chance constraints % % %
    
    for k = 0:N-1
        Fk = [zeros(nu,k*nu) eye(nu) zeros(nu,(N-k-1)*nu)];
        for i = 1:Nc
            ICDFu = norminv(1-deltau(i,k+1));
            controlConstraint = a_u(:,i)' * Fk * V + ICDFu * norm(S' * K' * Fk' * a_u(:,i)) <= b_u(i);
            Constraints = [Constraints,controlConstraint];
        end
    end
    
    % % % State chance constraints % % %
    
    Ac = PS.A_cone(1:3,1:3);
    cc = PS.c(1:3);
    H = [1 0 0; 0 0 1];
    Ip = [eye(3), zeros(3)];
    d = PS.d;
    for k = 1:N
        Ek = [zeros(nx,k*nx) eye(nx) zeros(nx,(N-k)*nx)];
        EXk = Ek * EX;
        xi_bar_k = H * Ac * Ip * Ek * EX;
        Rk = (cc' * Ip * EXk + d) - norm(xi_bar_k);
        stateConstraint = sqrt(2 * log(1 / deltax(k))) * norm(S' * IplusBK' * Ek' * Ip' * Ac' * H') <= Rk;
        Constraints = [Constraints,stateConstraint];
    end
end