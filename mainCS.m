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

[SolverTime,Jstar,Vsol,Ksol] = SolveCovarianceSteering(PS);

disp(Vsol(1,1));

% Calculate true deltas
Ns = PS.Ns;
N = PS.N;
deltabar = zeros(Ns,N);
for k = 1:N
    for j = 1:Ns
        deltabarjk = computeTrueDelta(PS,Vsol,Ksol,j,k);
        deltabar(j,k) = deltabarjk;
    end
end

disp(deltabar)

figure(1); 
subplot(1,2,1); hold on; grid on;
plot(PS.dt * [1:N],PS.deltax(2,:),'.k','MarkerSize',20);
plot(PS.dt * [1:N],PS.deltax(3,:),'ob','MarkerSize',10,'LineWidth',1);
ylim([0 6E-04]);
legend({'$\delta_r$','$\delta_u$'},'Interpreter','latex','FontSize',15,'Location','SouthWest');
xlabel('Time (s)','Interpreter','latex','FontSize',15);
ylabel('Allocated Risk','FontSize',15);
set(gca,'FontSize',15);
subplot(1,2,2); hold on; grid on;
plot(PS.dt * [1:N],deltabar(2,:),'.k','MarkerSize',20);
plot(PS.dt * [1:N],deltabar(3,:),'ob','MarkerSize',10,'LineWidth',1);
legend({'$\bar{\delta}_r$','$\bar{\delta}_u$'},'Interpreter','latex','FontSize',15,'Location','NorthWest');
xlabel('Time (s)','Interpreter','latex','FontSize',15);
ylabel('True Risk','FontSize',12);
set(gca,'FontSize',15);

PostProcess(2,100,Vsol,Ksol,PS);
PostProcess3D(3,100,Vsol,Ksol,PS,1);
V_N = computeVolumeN(PS,Ksol);