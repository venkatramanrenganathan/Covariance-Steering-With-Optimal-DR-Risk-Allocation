%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file solves the OCS problem with state DR Risk constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (C) Venkatraman Renganathan, Automatic Control LTH, Lund University 2022.
% (C) Joshua Pilipovsky, DCSL, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [time,obj,Vsol,Ksol] = SolveDRCovarianceSteering(PS, dynamicsSelectFlag)   
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
    % DR Risk Constraints
    Constraints = DRRiskConstraints(Constraints,EX,S,IplusBK,V,K,PS, dynamicsSelectFlag);
    
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

function Constraints = DRRiskConstraints(Constraints,EX,S,IplusBK,V,K,PS, dynamicsSelectFlag)
    a_x = PS.alphax;
    a_u = PS.alphau;
    b_x = PS.betax;
    b_u = PS.betau;
    Ns = PS.Ns;
    Nc = PS.Nc;
    N = PS.N;
    deltax = PS.deltax;
    deltau = PS.deltau;
    nx = PS.nx;
    nu = PS.nu;
    
    % State chance constraints
    for k = 1:N
        Ek = [zeros(nx,k*nx) eye(nx) zeros(nx,(N-k)*nx)];
        for j = 1:Ns            
            DR_Scale_X = sqrt((1 - deltax(j,k))/(deltax(j,k)));
            stateConstraint = a_x(:,j)'*Ek*EX + DR_Scale_X*norm(S'*IplusBK'*Ek'*a_x(:,j)) <= b_x(j);
            Constraints = [Constraints,stateConstraint];
        end
    end
    
    % Control chance constraints
    if (dynamicsSelectFlag == 1)
        for k = 0:N-1
            Fk = [zeros(nu,k*nu) eye(nu) zeros(nu,(N-k-1)*nu)];
            for i = 1:Nc
                DR_Scale_U = sqrt((1 - deltau(i,k+1))/(deltau(i,k+1)));
                controlConstraint = a_u(:,i)'*Fk*V + DR_Scale_U*norm(S'*K'*Fk'*a_u(:,i)) <= b_u(i);
                Constraints = [Constraints,controlConstraint];
            end
        end
    end
end