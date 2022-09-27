%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file solves the covariance-steering problem using iterative risk 
% allocation with distributionally robust convex cone state risk constraint. 
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
    PS = loadProblemSetup3D_CC(N);
elseif (dynamicsSelectFlag == 2)
    % Not implemented till now!!!
    N = 20; 
    PS = loadProblemSetupDoubleIntegrator(N);
end

%% Solve problem with Reverse Union Bound approximation to Cone Constraints
eps = 1E-05;
tstart = tic;
[PS,Vstar,Kstar,fstar,JstarVec,niter] = IterativeRiskAllocation_ConeConstraints_RUB(PS,eps, dynamicsSelectFlag, riskSelectFlag);
t_sim = toc(tstart);

%% Plot optimal cost for each IRA iteration
fignum = 1;
figure(fignum);
grid on;
plot(JstarVec,'.--k','MarkerSize',15);
xlabel('Iteration Number','FontSize',12);
ylabel('Optimal Cost','FontSize',12);
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 40);

% Plot The Post Processed Results
deltax_RUB_IRA = PS.deltax';
deltaTrue_RUB_IRA = computeTrueDeltaC_ConeConstraints_RUB(PS,Vstar,Kstar,fstar, riskSelectFlag);  
deltaTrue_RUB_IRA_max = max(max(deltaTrue_RUB_IRA));
deltaTrue_RUB_IRA_max = deltaTrue_RUB_IRA_max(:);
fignum = fignum + 1;
figure(fignum);
hold on;
plot(deltax_RUB_IRA - deltaTrue_RUB_IRA_max,'k','LineWidth',2);
xlabel('Iteration Number');
ylabel('$\bar{\delta}_{i,k} - \max_{i} \bar{\delta}_{i,k}$', 'interpreter', 'latex');
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 40);
hold off;

% Perform Monte-Carlo in Post Processing
if (dynamicsSelectFlag == 1)
    fignum = fignum + 1;
    PostProcess_ConeConstraints_New(fignum,100,Vstar,Kstar,PS)
    fignum = fignum + 1;
    PostProcess3D_CC(fignum,100,Vstar,Kstar,PS,1);
elseif(dynamicsSelectFlag == 2)
    fignum = fignum + 1;
    PostProcessDoubleIntegrator(fignum,100,Vstar,Kstar,PS);
end