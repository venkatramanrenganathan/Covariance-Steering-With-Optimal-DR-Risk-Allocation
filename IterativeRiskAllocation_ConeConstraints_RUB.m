function [PS,Vstar,Kstar,fstar,JstarVec,niter] = IterativeRiskAllocation_ConeConstraints_RUB(PS,eps)
% Risk Allocation / Covariance Steering Algorithm
niter = 0;
Jprev = Inf;
Jstar = 0;
JstarVec = [];

while abs(Jprev - Jstar) > eps
    alpha = 0.7 * 0.98 ^ niter; % Constraint tightening constant
    Jprev = Jstar;

    % Solve Lower Stage Covariance Steering %
    [SolverTime,Jstar,Vstar,Kstar,fstar] = SolveCovarianceSteering_ConeConstraints_RUB(PS);

    % Check active or inactive constraints
    [bool,diff,deltaTrue_RUB_max] = isActive_ConeConstraints_RUB(Vstar,Kstar,fstar,PS);

    %{
    % Look at distribution of constraints after initial "uniform" allocation
    % figure; hold on; grid on;
    % plot(diff(1,:),'.k','MarkerSize',15);
    % plot(diff(2,:),'.b','MarkerSize',15);
    % legend('Lower Constraint','Upper Constraint');
    %}

    Nactive = sum(sum(bool));
    if Nactive == 0 || Nactive == PS.N
        break
    end

    % Tighten INACTIVE constraints
    PS = tightenConstraints_ConeConstraints(PS,alpha,bool,deltaTrue_RUB_max);
    
    % Calculate residual probability (safety) margin
    deltaRes = PS.Deltax - sum(PS.deltax);
    
    % Loosen ACTIVE constraints
    PS = loosenConstraints_ConeConstraints(PS,Nactive,deltaRes,bool);

    % Iterate and store cost
    niter = niter + 1;
    JstarVec = [JstarVec,Jstar];
end


end