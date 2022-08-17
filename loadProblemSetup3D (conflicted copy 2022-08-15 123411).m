%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the problem setup (PS) given the length of the
% horizon N.
% (c) 2020 Joshua Pilipovsky, DCSL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PS = loadProblemSetup3D(N)
    PS.N = N;
    % Initial condition % [-1.15;-1.25;1.1;0;0;0];
    %PS.mu0 =  scale * [0.75;-1.0;0.75;0;0;0];
    PS.mu0 = [90; -120; 90; 0; 0; 0];
    %PS.mu0 = [100; -120; 90; 0; 0; 0];
    %PS.mu0 = [0.5;-1.5;0.75;0;0;0];
    PS.Sigma0 = 0.05*blkdiag(10,10,10,1,1,1);
    % Terminal Condition
    PS.muf = [0;0;0;0;0;0];
    PS.Sigmaf = PS.Sigma0; % 0.25 * PS.Sigma0;
    % Linear System
    dt = 4; 
    PS.dt = dt;
    PS = defineDynamics(PS, dt);
    % Probability threshold 
    PS.Ns = 4;
    PS.alphax = [-1 1 0 0 0 0;
                  0 1 1 0 0 0;
                  1 1 0 0 0 0;
                  0 1 -1 0 0 0]';
    PS.stateCC_offset = 20; % 6; % meters
    PS.betax = PS.stateCC_offset * ones(4,1);
    failProb = 0.1;
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
    PS.Deltau = 0.15; 
    % Initialize with uniform risk allocation
    PS = assignRisk(PS);
    % Objective function  
    PS = makeCost(PS, blkdiag(10,10,10,1,1,1),10^3*blkdiag(1,1,1));
    %PS = makeCost(PS, 10*blkdiag(10,10,10,1,1,1),10^4*blkdiag(1,1,1));
    % Make Large Matrices
    PS = makeLargeMatrices(PS);
end

function PS = defineDynamics(PS, dt)
    mu = 3.986E05; % [km]^3 / [s]^2
    Re = 6378.1; % [km]
    R0 = Re + 800; % Orbital radius for LEO: [km]
    omega = sqrt(mu/R0^3);
    mc = 100;
    PS.nx = 6;
    PS.nu = 3;
    PS.nw = 6;    
    [Ad,Bd] = get6dCwhStateAndInputMatrices(dt,omega,mc);
    PS.A = Ad;
    PS.B = Bd;
    PS.D = blkdiag(1E-04,1E-04,1E-04,5E-08,5E-08,5E-08);
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

