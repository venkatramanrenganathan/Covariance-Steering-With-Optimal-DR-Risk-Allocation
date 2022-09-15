%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file solves the covariance-steering problem under distributionally 
% robust state risk constraints. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (C) Venkatraman Renganathan, Automatic Control LTH, Lund University 2022.
% (C) Joshua Pilipovsky, DCSL, Georgia Tech 2022,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make a fresh start
close all; clear all; clc;

% Select the required flags
riskSelectFlag = 1; % 1 for Gaussian chance constraint, 2 for DR Risk Constraint
dynamicsSelectFlag = 1; % 1 for 3D spacecraft, 2 for Double Integrator 

% Set Horizon and load the data
if (dynamicsSelectFlag == 1)
    N = 15; 
    PS = loadProblemSetup3D(N);
elseif (dynamicsSelectFlag == 2)
    N = 20; 
    PS = loadProblemSetupDoubleIntegrator(N);
end

% Solve the Covariance Steering with DR Risk Constraints
if(riskSelectFlag == 1)
    [SolverTime,Jstar,Vsol,Ksol] = SolveCovarianceSteering(PS, dynamicsSelectFlag);
elseif(riskSelectFlag == 2)    
    [SolverTime,Jstar,Vsol,Ksol] = SolveDRCovarianceSteering(PS, dynamicsSelectFlag);
end

% Calculate true deltas
Ns = PS.Ns;
N = PS.N;
deltabar = zeros(Ns,N);
for k = 1:N
    for j = 1:Ns
        if(riskSelectFlag == 1)
            deltabarjk = computeTrueDelta(PS,Vsol,Ksol,j,k);
        elseif(riskSelectFlag == 2)
            deltabarjk = computeDRTrueDelta(PS,Vsol,Ksol,j,k);
        end
        deltabar(j,k) = deltabarjk;
    end
end

%% Plotting 
figNum = 1;
figure(figNum); 
subplot(1,2,1); 
hold on; 
grid on;
plot(PS.dt*[1:N],PS.deltax(2,:),'*k','MarkerSize',20);
plot(PS.dt*[1:N],PS.deltax(3,:),'or','MarkerSize',10,'LineWidth',1);
ylim([0 9E-04]);
legend({'$\delta_r$','$\delta_u$'},'Interpreter','latex','FontSize',30,'Location','SouthWest');
xlabel('Time (s)','Interpreter','latex','FontSize',15);
ylabel('Allocated Risk','Interpreter','latex','FontSize',15);
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 2);
set(a, 'linewidth', 2);
set(a, 'FontSize', 40);
subplot(1,2,2); 
hold on; 
grid on;
plot(PS.dt*[1:N],deltabar(2,:),'*k','MarkerSize',20);
plot(PS.dt*[1:N],deltabar(3,:),'or','MarkerSize',10,'LineWidth',1);
ylim([0 9E-04]);
legend({'$\bar{\delta}_r$','$\bar{\delta}_u$'},'Interpreter','latex','FontSize',30,'Location','SouthWest');
xlabel('Time (s)','Interpreter','latex','FontSize',15);
ylabel('True Risk','Interpreter','latex','FontSize',12);
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 2);
set(a, 'linewidth', 2);
set(a, 'FontSize', 40);

if (dynamicsSelectFlag == 1)
    figNum = figNum + 1;
    PostProcess(figNum,100,Vsol,Ksol,PS);
    figNum = figNum + 1;
    PostProcess3D(figNum,100,Vsol,Ksol,PS,1);
    V_N = computeVolumeN(PS,Ksol);
elseif(dynamicsSelectFlag == 2)
    figNum = figNum + 1;
    PostProcessDoubleIntegrator(figNum,100,Vsol,Ksol,PS);
end