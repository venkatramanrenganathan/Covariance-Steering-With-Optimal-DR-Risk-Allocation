%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the problem setup (PS) given the length of the
% horizon N.
% (c) 2020 Joshua Pilipovsky, DCSL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PS = loadProblemSetup(N)
    PS.N = N;
    % Initial condition
    PS.mu0 =  [-1.15;-1.25;0;0];
    PS.Sigma0 = 1E-02 * blkdiag(0.1,0.1,0.01,0.01);
    % Terminal Condition
    PS.muf = [0;0;0;0];
    PS.Sigmaf = 0.25 * PS.Sigma0;
    % Linear System
    dt = 1;
    PS = defineDynamics(PS, dt);
    % Probability threshold 
    PS.M = 2;
    PS.alpha = [-1 1 0 0;
                1 1 0 0]';
    PS.beta = [0.05; 0.05];
    PS.Delta = 0.02; % JOINT probability of failure 
    % Initialize with uniform risk allocation
    PS = assignRisk(PS);
    % Objective function  
    PS = makeCost(PS, blkdiag(10,10,1,1),blkdiag(10^3,10^3));
    % Make Large Matrices
    PS = makeLargeMatrices(PS);
end

function PS = defineDynamics(PS, dt)
    mu = 3.986E05; % [km]^3 / [s]^2
    Re = 6378.1; % [km]
    R0 = Re + 850; % Orbital radius for LEO: [km]
    omega = sqrt(mu/R0^3);
    PS.nx = 4;
    PS.nu = 2;
    PS.nw = 4;    
%     Ac = [0 1 0 0;
%             3*omega^2 0 2*omega 0;
%             0 0 0 1;
%             0 -2*omega 0 0];
    Ac = [0 0 1 0;
          0 0 0 1;
          3*omega^2 2*omega 0 0;
          0 0 -2*omega 0];
    PS.A = eye(PS.nx) + dt * Ac; % First order approx.
    Bc =  (1/300) * [0 0;
          0 0;
          1 0;
          0 1];
%     Bc = [0 0;
%           1 0;
%           0 0;
%           0 1];
    PS.B = dt * Bc; % First order approx.
    PS.D = blkdiag(1E-04,1E-04,5E-08,5E-08);
end

function PS = assignRisk(PS)
    N = PS.N;
    M = PS.M;
    Delta = PS.Delta;
    delta = zeros(M,N);
    for k = 1:N
        for j = 1:M
            delta(j,k) = Delta / (N * M);
        end
    end
    PS.delta = delta;
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

