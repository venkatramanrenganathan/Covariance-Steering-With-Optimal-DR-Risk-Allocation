% This script compares the scaling constant for distributionally robust
% risk constraint and the Gaussian chance constraint.

% make a fresh start
close all; clear all; clc;

% make a vector of risk parameters
deltas = 0.001:0.001:0.50;

% Infer the length of deltas vector
numDeltas = size(deltas,2);

% Placeholder to store scalings
GaussianScalings = zeros(numDeltas, 1);
GaussianTailings = zeros(numDeltas, 1);
distrRobScalings = zeros(numDeltas, 1);
distrRobTailings = zeros(numDeltas, 1);

% Loop and populate the scalings
for i = 1:numDeltas
    % GaussianScalings(i,1) = sqrt(2)*erfinv(1-2*deltas(i));
    GaussianScalings(i,1) = norminv(1-deltas(i), 0, 1);
    GaussianTailings(i,1) = norminv(deltas(i), 0, 1);
    distrRobScalings(i,1) = sqrt((1-deltas(i))/(deltas(i)));
    distrRobTailings(i,1) = sqrt((deltas(i))/(1-deltas(i)));
end

% Plot the results
figure(1)
plot(deltas, GaussianScalings, '-r');
hold on;
plot(deltas, GaussianTailings, '-.r');
hold on;
plot(deltas, distrRobScalings, '-b');
hold on;
plot(deltas, distrRobTailings, '-.b');
legend('$\Phi^{-1}(1-\delta)$', '$\Phi^{-1}(\delta)$', '$\mathcal{Q}(\delta) = \sqrt{\frac{1-\delta}{\delta}}$', '$\mathcal{Q}^{-1}(\delta) = \sqrt{\frac{\delta}{1-\delta}}$', 'interpreter', 'latex');
xlabel('Risk parameter $\delta$', 'interpreter', 'latex');
ylabel('Scaling', 'interpreter', 'latex');
hold off;
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 40);