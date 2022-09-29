%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file solves the covariance-steering problem using iterative risk 
% allocation with distributionally robust polytopic state risk constraints. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (C) Venkatraman Renganathan, Automatic Control LTH, Lund University 2022,
% (C) Joshua Pilipovsky, DCSL, Georgia Tech 2022,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make a fresh start
close all; clear all; clc;

% Select the required flags
riskSelectFlag = 2; % 1 for Gaussian chance constraint, 2 for DR Risk Constraint
dynamicsSelectFlag = 1; % 1 for 3D spacecraft, 2 for Double Integrator 

% Set Horizon and load the data
if (dynamicsSelectFlag == 1)
    N = 15; 
    eps = 1E02;
    PS = loadProblemSetup3D(N);
elseif (dynamicsSelectFlag == 2)
    N = 20; 
    eps = 1E02;
    PS = loadProblemSetupDoubleIntegrator(N);
end

% Run iterative risk allocation algorithm
[PS,Vstar,Kstar,JstarVec,iter] = IterativeRiskAllocation(PS,eps, dynamicsSelectFlag, riskSelectFlag);

% Plot optimal cost history
figNum = 1;
figure(figNum); 
grid on;
plot(JstarVec,'.--k','MarkerSize',15);
xlim([0 iter]);
xlabel('Iteration Number');
ylabel('Optimal Cost');
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 6);
set(a, 'linewidth', 6);
set(a, 'FontSize', 60);

% Calculate true deltas after IRA
Ns = PS.Ns;
N = PS.N;
PS.deltabar = zeros(Ns,N);
for k = 1:N
    for j = 1:Ns
        if(riskSelectFlag == 1)
            % Gaussian case
            PS.deltabar(j,k) = computeTrueDelta(PS,Vstar,Kstar,j,k);
        elseif(riskSelectFlag == 2)
            % DR case
            PS.deltabar(j,k) = computeDRTrueDelta(PS,Vstar,Kstar,j,k);
        end
    end
end

% Compute Volume and Fuel Cost for RPO Problem
V_N = computeVolumeN(PS,Kstar);
FuelCosts = norm(Vstar,1);

%%  Plot results for RPO Problem
figNum = figNum + 1;
figure(figNum); 
subplot(1,2,1); hold on; grid on;
plot(PS.dt * [1:PS.N],PS.deltax(2,:),'*k','MarkerSize',30);
if (dynamicsSelectFlag == 1)
    plot(PS.dt * [1:PS.N],PS.deltax(3,:),'or','MarkerSize',20,'LineWidth',1);
    legend({'$\delta_r$','$\delta_u$'},'Interpreter','latex','FontSize',20);
elseif(dynamicsSelectFlag == 2)
    legend({'$\delta_r$'},'Interpreter','latex','FontSize',20);
end
xlabel('Time (s)','Interpreter','latex','FontSize',15);
ylabel('Allocated Risk','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
subplot(1,2,2); hold on; grid on;
plot(PS.dt * [1:PS.N],PS.deltabar(2,:),'*k','MarkerSize',30);
if (dynamicsSelectFlag == 1)
    plot(PS.dt * [1:PS.N],PS.deltabar(3,:),'or','MarkerSize',20,'LineWidth',1);
    legend({'$\bar{\delta}_r$','$\bar{\delta}_u$'},'Interpreter','latex','FontSize',20);
elseif(dynamicsSelectFlag == 2)
    legend({'$\bar{\delta}_r$'},'Interpreter','latex','FontSize',20);
end
xlabel('Time (s)','Interpreter','latex','FontSize',15);
ylabel('True Risk','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 60);

% Perform Monte-Carlo in Post Processing
numMonteTrials = 500;
if (dynamicsSelectFlag == 1)
    figNum = figNum + 1;
    PostProcess(figNum,numMonteTrials,Vstar,Kstar,PS);
%     figNum = figNum + 1;
%     PostProcess3D(figNum,numMonteTrials,Vstar,Kstar,PS,1);
%     figNum = figNum + 2; % Because PostProcess3D produces 2 figures
%     PostProcess3D(7,numMonteTrials,Vstar,Kstar,PS,0);
elseif (dynamicsSelectFlag == 2)
    figNum = figNum + 1;
    PostProcessDoubleIntegrator(figNum,numMonteTrials,Vstar,Kstar,PS);
end
