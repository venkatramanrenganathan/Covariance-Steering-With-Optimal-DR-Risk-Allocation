function [PS,Vsol,Ksol,JstarVec,niter] = IterativeRiskAllocation(PS,eps, dynamicsSelectFlag, riskSelectFlag)
    % Risk Allocation / Covariance Steering Algorithm
    niter = 0;
    Jprev = Inf;
    Jstar = 0;
    JstarVec = [];
    
    %%% BEGIN IRA LOOP %%%
    while (abs(Jprev - Jstar) > eps)
        alpha = 0.7 * 0.98 ^ niter; % Constraint tightening constant
        Jprev = Jstar;
        
        fprintf('Risk Allocation Iteration No: %d \n', niter+1);

        % Solve Lower Stage Covariance Steering with DR Risk/ Chance Constraints
        if(riskSelectFlag == 1)
            [SolverTime,Jstar,Vsol,Ksol] = SolveCovarianceSteering(PS, dynamicsSelectFlag);
        elseif(riskSelectFlag == 2)    
            [SolverTime,Jstar,Vsol,Ksol] = SolveDRCovarianceSteering(PS, dynamicsSelectFlag);
        end
        
        % Calculate true deltas after Uniform Allocation (Baseline)
        Ns = PS.Ns;
        N = PS.N;
        PS.deltabar0 = zeros(Ns,N);
        for k = 1:N
            for j = 1:Ns
                if(riskSelectFlag == 1)
                    PS.deltabar0(j,k) = computeTrueDelta(PS,Vsol,Ksol,j,k);
                elseif(riskSelectFlag == 2)
                    PS.deltabar0(j,k) = computeDRTrueDelta(PS,Vsol,Ksol,j,k);
                end
            end
        end

        % Check active or inactive constraints
        [bool,diff] = isActive(Vsol,Ksol,PS, riskSelectFlag);        

        Nactive = sum(sum(bool));
        if Nactive == 0 || Nactive == PS.N * PS.Ns
            break
        end

        % Tighten INACTIVE constraints
        PS = tightenConstraints(PS,alpha,Vsol,Ksol,bool, riskSelectFlag);
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