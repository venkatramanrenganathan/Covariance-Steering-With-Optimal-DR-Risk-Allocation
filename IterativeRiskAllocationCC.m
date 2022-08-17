function [PS,Vsol,Ksol,JstarVec,niter] = IterativeRiskAllocationCC(PS,eps,char)
    % Risk Allocation / Covariance Steering Algorithm
    niter = 0;
    Jprev = Inf;
    Jstar = 0;
    JstarVec = [];
    
    if strcmp(char,'GEO')
    %%% BEGIN IRA LOOP %%%
        while abs(Jprev - Jstar) > eps
            alpha = 0.7 * 0.98 ^ niter; % Constraint tightening constant
            Jprev = Jstar;

            % Solve Lower Stage Covariance Steering %
            [SolverTime,Jstar,Vsol,Ksol] = SolveCovarianceSteeringCC(PS,'GEO');

            % Check active or inactive constraints
            [bool,diff] = isActive_CC(Vsol,Ksol,PS,'GEO');

            %{
            % Look at distribution of constraints after initial "uniform" allocation
            % figure; hold on; grid on;
            % plot(diff(1,:),'.k','MarkerSize',15);
            % plot(diff(2,:),'.b','MarkerSize',15);
            % legend('Lower Constraint','Upper Constraint');
            %}

            Nactive = sum(bool);
            if Nactive == 0 || Nactive == PS.N
                break
            end

            % Tighten INACTIVE constraints
            PS = tightenConstraints_CC(PS,alpha,Vsol,Ksol,bool,'GEO');
            % Calculate residual probability (safety) margin
            deltaRes = PS.Deltax - sum(PS.deltax);
            % Loosen ACTIVE constraints
            PS = loosenConstraints_CC(PS,Nactive,deltaRes,bool,'GEO');

            % Iterate and store cost
            niter = niter + 1;
            JstarVec = [JstarVec,Jstar];
        end
    %%% END IRA LOOP %%%
    else
        while abs(Jprev - Jstar) > eps
            alpha = 0.7 * 0.98 ^ niter; % Constraint tightening constant
            Jprev = Jstar;

            % Solve Lower Stage Covariance Steering %
            [SolverTime,Jstar,Vsol,Ksol] = SolveCovarianceSteeringCC(PS,'TS');

            % Check active or inactive constraints
            [bool,diff] = isActive_CC(Vsol,Ksol,PS,'TS');

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
            PS = tightenConstraints_CC(PS,alpha,Vsol,Ksol,bool,'TS');
            % Calculate residual probability (safety) margin
            deltaRes = PS.Deltax - sum(PS.deltax);
            % Loosen ACTIVE constraints
            PS = loosenConstraints_CC(PS,Nactive,deltaRes,bool,'TS');

            % Iterate and store cost
            niter = niter + 1;
            JstarVec = [JstarVec,Jstar];
        end
    end
        
end