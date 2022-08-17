%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file solves the covariance-steering problem under state chance 
% constraints. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Joshua Pilipovsky, DCSL, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping
close all;
clear all;
clc;

%% Setup Problem 
% Set Horizon
N = 15;
PS = loadProblemSetup3D_CC(N);

% Plot uncontrolled trajectories
%{
MCnum = 50;
nx = PS.nx;
nu = PS.nu;
nw = PS.nw;
N = PS.N;
figure; 
for mc = 1:MCnum
   x0_MC = PS.mu0 + sqrtm(PS.Sigma0)*randn(nx,1);
   x_MC = zeros(nx,N+1);
   x_MC(:,1) = x0_MC;
   for k = 1:N
       w = randn(nw,1);
       x_MC(:,k+1) = PS.A*x_MC(:,k) + PS.D*w;
   end       
   plot3(x_MC(1,:),x_MC(2,:),x_MC(3,:),'color',[0.5,0.5,0.5]);
   hold on; grid on;
end
hold on; xlabel('x'); ylabel('y'); zlabel('z');
plot3(0,0,0,'.r','MarkerSize',10);
%}

%% Solve problem with Reverse Union Bound approximation to Cone Constraints
eps = 1E-05;
tstart = tic;
[PS,Vstar,Kstar,fstar,JstarVec,niter] = IterativeRiskAllocation_ConeConstraints_RUB(PS,eps);
t_sim = toc(tstart);
%[SolverTime,Jstar_RUB,Vstar,Kstar,fstar] = SolveCovarianceSteering_ConeConstraints_RUB(PS);


 %% Plot optimal cost for each IRA iteration
% Plot optimal cost history
figure; grid on;
plot(JstarVec,'.--k','MarkerSize',15);
xlabel('Iteration Number','FontSize',12);
ylabel('Optimal Cost','FontSize',12);

%% Post Process Results

deltax_RUB_IRA = PS.deltax';
deltaTrue_RUB_IRA = computeTrueDeltaC_ConeConstraints_RUB(PS,Vstar,Kstar,fstar);  
deltaTrue_RUB_IRA_max = max(max(deltaTrue_RUB_IRA));
deltaTrue_RUB_IRA_max = deltaTrue_RUB_IRA_max(:);
figure; hold on;
plot(deltax_RUB_IRA - deltaTrue_RUB_IRA_max,'k','LineWidth',2);

PostProcess_ConeConstraints_New(1,100,Vstar,Kstar,PS)
PostProcess3D_CC(3,100,Vstar,Kstar,PS,1);

%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Solve problem with Three-Cut approximation to Cone Constraints
eps = 1E-05;
[PS,Vstar,Kstar,fstar,tstar,JstarVec,niter] = IterativeRiskAllocation_ConeConstraints_TC(PS,eps);
%[SolverTime,Jstar_TC,Vstar,Kstar,fstar,tstar] = SolveCovarianceSteering_ConeConstraints_ThreeCut(PS);

%% Plot optimal cost for each IRA iteration
% Plot optimal cost history
figure; grid on;
plot(JstarVec,'.--k','MarkerSize',15);
xlabel('Iteration Number','FontSize',12);
ylabel('Optimal Cost','FontSize',12);

%% Post Process Results 
%deltaTrue_TC = computeTrueDeltaC_ConeConstraints_ThreeCut(PS,Vstar,fstar,tstar);

deltax_TC_IRA = PS.deltax';
deltaTrue_TC_IRA = computeTrueDeltaC_ConeConstraints_ThreeCut(PS,Vstar,fstar,tstar);  
deltaTrue_TC_IRA_max = max(max(deltaTrue_TC_IRA));
deltaTrue_TC_IRA_max = deltaTrue_TC_IRA_max(:);
figure; hold on;
plot(deltax_TC_IRA - deltaTrue_TC_IRA_max,'k','LineWidth',2);

PostProcess_ConeConstraints_New(1,100,Vstar,Kstar,PS)
PostProcess3D_CC(3,100,Vstar,Kstar,PS,1);

%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Solve Problem with Geometric Cone Constraints
% Run iterative risk allocation algorithm
eps = 10;
[PS,Vstar,Kstar,JstarVec,iter] = IterativeRiskAllocationCC(PS,eps,'GEO');
%[SolverTime,Jstar_GEO,Vstar,Kstar] = SolveCovarianceSteeringCC(PS,'GEO_new');

PostProcess_ConeConstraints_New(1,100,Vstar,Kstar,PS)
PostProcess3D_CC(3,100,Vstar,Kstar,PS,1);

 %% Plot optimal cost for each IRA iteration
% Plot optimal cost history
figure; grid on;
plot(JstarVec,'.--k','MarkerSize',15);
xlabel('Iteration Number','FontSize',12);
ylabel('Optimal Cost','FontSize',12);

%% Calculate True Risks
PS.deltabar_GEO = zeros(1,N);
for k = 1:N
    deltabark = computeTrueDeltaCC2(PS,Vstar,Kstar,k);
    PS.deltabar_GEO(k) = deltabark;
end

%% Plot Results
PostProcessCC(2,100,Vstar,Kstar,PS);
PostProcess3D_CC(3,100,Vstar,Kstar,PS,1);
% Final Volume
VN_GEO = computeVolumeN(PS,Kstar);
% Total Fuel Costs
FuelCosts_GEO = norm(Vstar,1);

%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Setup Problem 
clear all; clc;
% Set Horizon
N = 15;
PS = loadProblemSetup3D_CC(N);
%% Solve Problem with Two-Sided Cone Constraints
% Run iterative risk allocation algorithm
eps = 1E02;
[PS,Vstar,Kstar,JstarVec,iter] = IterativeRiskAllocationCC(PS,eps,'TS');
%[SolverTime,Jstar_TS,Vsol,Ksol] = SolveCovarianceSteeringCC(PS,'TS');

%% Calculate True Risks
PS.deltabar_TS = zeros(4,N);
for k = 1:N
    deltabark = computeTrueDeltaCC(PS,Vstar,Kstar,k);
    PS.deltabar_TS(:,k) = deltabark;
end

%% Plot Results
PostProcessCC(6,100,Vstar,Kstar,PS);
PostProcess3D_CC(7,100,Vstar,Kstar,PS,1);
% Final Volume
VN_TS = computeVolumeN(PS,Kstar);
% Total Fuel Costs
FuelCosts_TS = norm(Vstar,1);

%% Plot optimal cost for each IRA iteration
% Plot optimal cost history
figure; grid on;
plot(JstarVec,'.--k','MarkerSize',15);
xlabel('Iteration Number','FontSize',12);
ylabel('Optimal Cost','FontSize',12);


