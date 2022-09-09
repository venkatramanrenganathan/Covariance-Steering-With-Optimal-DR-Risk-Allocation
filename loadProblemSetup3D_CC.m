%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the problem setup (PS) given the length of the
% horizon N.
% (c) 2022 Joshua Pilipovsky, DCSL
% (c) 2022 Venkatraman Renganathan, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PS = loadProblemSetup3D_CC(N)
    
    % Time horizon
    PS.N = N;
    
    % Initial Conditions
    PS.mu0 = [10; 120; 90; 0; 0; 0];
    PS.Sigma0 = 0.01*blkdiag(10,10,10,1,1,1);
    
    % Infer the system size
    n = size(PS.mu0,1);
    
    % Terminal Conditions
    PS.muf = zeros(n,1);
    PS.Sigmaf = 0.25*PS.Sigma0;
    
    % Linear System dynamics
    dt = 4;
    PS.dt = dt;
    PS = defineDynamics(PS, dt);
    
    % Cone state chance constraints
    PS.Deltax = 0.05; % 0.03; % JOINT probability of failure 
    
    % Cone parametrization: || A * x + b || <= c'*xk + d    
    psi          = -atand(3/4);
    coneAngle    = 15;
    gamma        = tand(coneAngle);
    PS.psi       = psi;
    PS.gamma     = gamma;
    PS.coneAngle = coneAngle;
    
    
    Rx = [1 0         0; 
          0 cosd(psi) -sind(psi); 
          0 sind(psi) cosd(psi)];
    Rxbar    = blkdiag(Rx,Rx);
    PS.Rxbar = Rxbar;
    
    A_x       = blkdiag(1,0,1,0,0,0);
    A_cone_xy = [cosd(psi)^2 -sind(psi) * cosd(psi);
                 -sind(psi) * cosd(psi) sind(psi)^2];
    A_cone_xyz = blkdiag(A_cone_xy,1);
    PS.A_cone  = A_x * Rxbar;
    PS.b = zeros(n,1);
    PS.c = Rxbar'*[0; gamma; zeros(n-2,1)];
    PS.d = 10;
    PS.stateCC_offset = PS.d; % m
    
    % Two sided chance constraint parameter
    PS.beta = 1/PS.nx;
    
    % Polyhedron control chance constraints
    PS.Nc = 6;
    PS.alphau = [1 0 0;
                 0 1 0;
                 0 0 1;
                 -1 0 0;
                 0 -1 0;
                 0 0 -1]';
    PS.umax = 30; % Newtons
    PS.betau = PS.umax * ones(PS.Nc,1);
    PS.Deltau = 0.10; 
    
    % Initialize with uniform risk allocation
    PS = assignRisk(PS);
    % Objective function  
    PS = makeCost(PS, blkdiag(10,10,10,1,1,1),blkdiag(10^3,10^3,10^3));
    % Make Large Matrices
    PS = makeLargeMatrices(PS);
end

function PS = defineDynamics(PS, dt)
    mu = 3.986E05; % [km]^3 / [s]^2
    Re = 6378.1; % [km]
    R0 = Re + 800; % Orbital radius for LEO: [km]
    omega = sqrt(mu/R0^3);
    mc    = 100; % [kg]
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
    Nc = PS.Nc;
    Deltax = PS.Deltax;
    Deltau = PS.Deltau;
    deltax = zeros(1,N);
    deltau = zeros(Nc,N);
    for k = 1:N
        deltax(k) = Deltax / N;
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

