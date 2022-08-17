%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file solves the covariance-steering problem under state chance 
% constraints. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Joshua Pilipovsky, DCSL, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;

% Set Horizon
N = 15;
PS = loadProblemSetup3D(N);

% Run iterative risk allocation algorithm
eps = 1E02;
[PS,Vstar,Kstar,JstarVec,iter] = IterativeRiskAllocation(PS,eps);

% Plot optimal cost history
figure(1); grid on;
plot(JstarVec,'.--k','MarkerSize',15);
xlabel('Iteration Number','FontSize',12);
ylabel('Optimal Cost','FontSize',12);

% Calculate true deltas after IRA
Ns = PS.Ns;
N = PS.N;
PS.deltabar = zeros(Ns,N);
for k = 1:N
    for j = 1:Ns
        deltabarjk = computeTrueDelta(PS,Vstar,Kstar,j,k);
        PS.deltabar(j,k) = deltabarjk;
    end
end

% Plot results for RPO Problem
PlotResultsRPO(Vstar,Kstar,PS);
V_N = computeVolumeN(PS,Kstar);
FuelCosts = norm(Vstar,1);
