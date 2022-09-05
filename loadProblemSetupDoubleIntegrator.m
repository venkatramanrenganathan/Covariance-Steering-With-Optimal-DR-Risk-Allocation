%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the problem setup (PS) given the length of the
% horizon N for a double integrator system
% (C) 2022 Venkatraman Renganathan
%          Department of Automatic Control - LTH
%          Lund University, Sweden
%          Email: venkatraman.renganathan@control.lth.se
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PS = loadProblemSetupDoubleIntegrator(N)
    PS.N = N;
    % Initial condition 
    PS.mu0 = [-10; 1; 0; 0]; % [-10; 1; 0; 0];
    PS.Sigma0 = diag([0.1,0.1,0.01,0.01]); 
    % Terminal Condition
    PS.muf = [0;0;0;0];
    PS.Sigmaf = 0.25*PS.Sigma0; % 0.75*PS.Sigma0; 
    % Linear System
    dt = 0.2; 
    PS.dt = dt;
    PS = defineDynamics(PS, dt);
    % Probability threshold 
    PS.Ns = 2;
    PS.alphax = [0.2 -1 0 0;
                 0.2 1  0 0]';
    PS.alphax = PS.alphax;         
    PS.stateCC_offset = 0.2; % meters
    PS.betax = PS.stateCC_offset * ones(4,1);
    failProb = 0.10; % 0.0005;
    PS.Deltax = failProb; % 0.03; % JOINT probability of failure 
    % Control Chance Constraints
    PS.Nc = 6;
    PS.alphau = [1 0 0;
                 0 1 0;
                 0 0 1;
                 -1 0 0;
                 0 -1 0;
                 0 0 -1]';
    PS.umax = 30; % Newtons
    PS.betau = PS.umax * ones(PS.Nc,1);
    PS.Deltau = failProb; % 0.20; 
    
    % Two sided chance constraint parameter
    PS.beta = 1/size(PS.mu0, 1);
    
    % Initialize with uniform risk allocation
    PS = assignRisk(PS);
    % Objective function  
    PS = makeCost(PS, blkdiag(10,10,1,1),10^3*blkdiag(1,1));
    % Make Large Matrices
    PS = makeLargeMatrices(PS);
end

function PS = defineDynamics(PS, dt)
    dt2 = dt*dt;
    PS.A = [1 0 dt 0
            0 1 0  dt
            0 0 1  0
            0 0 0  1];
    PS.B = [dt2 0
            0   dt2
            dt  0
            0   dt];
    PS.nx = size(PS.A, 1);
    PS.nu = size(PS.B, 2);
    PS.nw = PS.nx;
    PS.D  = 0.001*eye(PS.nx);
end

function PS = assignRisk(PS)
    N = PS.N;
    Ns = PS.Ns;
    Nc = PS.Nc;
    Deltax = PS.Deltax;
    Deltau = PS.Deltau;
    deltax = zeros(Ns,N);
    deltau = zeros(Nc,N);
    for k = 1:N
        for j = 1:Ns
            deltax(j,k) = Deltax / (N * Ns);
        end
        for i = 1:Nc
            deltau(i,k) = Deltau / (N * Nc);
        end
    end
    PS.deltax = deltax;
    PS.deltau = deltau;
end

function PS = makeCost(PS,Q0,R0)
    Q = [];
    R = [];
    for i = 1:PS.N
        Q = blkdiag(Q,Q0);
        R = blkdiag(R,R0);
    end
    Q = blkdiag(Q,zeros(PS.nx));
    PS.Q = Q;
    PS.R = R;
end

function PS = makeLargeMatrices(PS)
    nx = PS.nx;
    nu = PS.nu;
    nw = PS.nw;
    N = PS.N;
    A = PS.A;
    B = PS.B;
    D = PS.D;
    ScriptA = eye(nx);
    ScriptB = zeros(nx,nu*N);
    ScriptD = zeros(nx,nw*N);
    for k = 1:N
        ScriptA = [ScriptA; A*ScriptA(end-nx+1:end,:)];
        ScriptB = [ScriptB; A*ScriptB(end-nx+1:end,1:(k-1)*nu) B zeros(nx,(N-k)*nu)];
        ScriptD = [ScriptD; A*ScriptD(end-nx+1:end,1:(k-1)*nw) D zeros(nx,(N-k)*nw)];
    end
    PS.ScriptA = ScriptA;
    PS.ScriptB = ScriptB;
    PS.ScriptD = ScriptD;
end

