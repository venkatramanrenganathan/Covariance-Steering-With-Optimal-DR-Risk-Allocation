function [PS,Vsol,Ksol,JstarVec,niter] = IterativeRiskAllocation(PS,eps)
    % Risk Allocation / Covariance Steering Algorithm
    niter = 0;
    Jprev = Inf;
    Jstar = 0;
    JstarVec = [];
    
    %%% BEGIN IRA LOOP %%%
    while abs(Jprev - Jstar) > eps
        alpha = 0.7 * 0.98 ^ niter; % Constraint tightening constant
        Jprev = Jstar;

        % Solve Lower Stage Covariance Steering %
        [SolverTime,Jstar,Vsol,Ksol] = SolveCovarianceSteering(PS);
        
        % Calculate true deltas after Uniform Allocation (Baseline)
        Ns = PS.Ns;
        N = PS.N;
        PS.deltabar0 = zeros(Ns,N);
        for k = 1:N
            for j = 1:Ns
                deltabarjk0 = computeTrueDelta(PS,Vsol,Ksol,j,k);
                PS.deltabar0(j,k) = deltabarjk0;
            end
        end

        % Check active or inactive constraints
        [bool,diff] = isActive(Vsol,Ksol,PS);

        % Look at distribution of constraints after initial "uniform" allocation
        % figure; hold on; grid on;
        % plot(diff(1,:),'.k','MarkerSize',15);
        % plot(diff(2,:),'.b','MarkerSize',15);
        % legend('Lower Constraint','Upper Constraint');

        Nactive = sum(sum(bool));
        if Nactive == 0 || Nactive == PS.N * PS.Ns
            break
        end

        % Tighten INACTIVE constraints
        PS = tightenConstraints(PS,alpha,Vsol,Ksol,bool);
        % Calculate residual probability (safety) margin
        deltaRes = PS.Deltax - sum(sum(PS.deltax));
        % Loosen ACTIVE constraints
        PS = loosenConstraints(PS,Nactive,deltaRes,bool);

        % Iterate and store cost
        niter = niter + 1;
        JstarVec = [JstarVec,Jstar];
    end
    %%% END IRA LOOP %%%
end