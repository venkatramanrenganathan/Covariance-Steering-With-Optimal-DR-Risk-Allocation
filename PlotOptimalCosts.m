%% Plot optimal costs over each IRA iteration
clear all; close all; clc;
str_TC = load('JstarVec_TC.mat');
str_RUB = load('JstarVec_RUB.mat');
str_GEO = load('JstarVec_GEO.mat');
str_poly = load('JstarVec_poly.mat');

JstarVec_TC = str_TC.JstarVec;
JstarVec_RUB = str_RUB.JstarVec;
JstarVec_GEO = str_GEO.JstarVec;
JstarVec_poly = str_poly.JstarVec;

figure; hold on; grid on;
plot(JstarVec_poly,'-k','LineWidth',2);
%plot(JstarVec_TC,'--k','LineWidth',2);
%plot(JstarVec_RUB,'-.k','LineWidth',2);
%plot(JstarVec_GEO,':k','LineWidth',2);
%legend({'Polyhedron','Three-Cut','RUB','Geometric'},'FontSize',15, ...,
%        'Location','NorthEast');
xlabel('Iteration','FontSize',15);
ylabel('Optimal Cost','FontSize',15);
set(gca,'FontSize',15);
