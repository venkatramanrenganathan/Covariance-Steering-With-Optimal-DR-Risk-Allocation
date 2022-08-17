function PostProcess_ConeConstraints_New(fignum,MCnum,V,K,PS)
    % Color for plotting MC's
    grey = [0.7 0.7 0.7];
    lightred = [255 154 162] / 255;
    
    % Plot cone cross section
    figure(fignum); hold on;
    %str1 = sprintf('%d * x.^2 > (%d * y + %d).^2',PS.eta,PS.gamma,PS.d);
    
    gamma = PS.gamma;
    psi = 0;
    d = PS.d;
    
% %     Xi = [gamma * cosd(psi) + sind(psi) cosd(psi) - gamma * sind(psi);
% %           gamma * cosd(psi) - sind(psi) cosd(psi) + gamma * sind(psi)];
% %       
% %     str1 = sprintf('%d * y < %d * x - %d',Xi(1,1),Xi(1,2),d);
% %     str2 = sprintf('%d * y < -%d * x - %d',Xi(2,1),Xi(2,2),d);
%     
%     str1 = sprintf('x.^2 + y.^2 * sin(%d)^2 > (%d * cos(%d) * y + %d)^2', ...
%                    psi, gamma, psi, d);
%     
%     ineqplot(str1,[-20 120 -20 160],lightred);
% %     ineqplot(str2,[-20 120 -20 160],lightred);
    
%     y_cone_outline = linspace(-20, 160, 1000);
%     %x_cone_outline = sqrt(PS.eta)^(-1) * (PS.gamma * y_cone_outline + PS.d);
%     x_cone_outline1 = Xi(1,2)^(-1) * (y_cone_outline * Xi(1,1) + d);
%     x_cone_outline2 = -Xi(2,2)^(-1) * (y_cone_outline * Xi(2,1) + d);
%     
%     plot(x_cone_outline1,y_cone_outline,'k','LineWidth',1);
%     plot(x_cone_outline2,y_cone_outline,'k','LineWidth',1);

    y = linspace(-20,160);
    f1_psi = gamma^2 * cosd(psi)^2 - sind(psi)^2;
    f2_psi = 2 * gamma * d * cosd(psi);
    x_proot = sqrt(f1_psi .* y.^2 + f2_psi .* y + d^2);
    x_mroot = -sqrt(f1_psi .* y.^2 + f2_psi .* y + d^2);
    
    plot(x_proot,y,'k','LineWidth',2);
    plot(x_mroot,y,'k','LineWidth',2);
    
    %xlim([-20 120]);
    ylim([-20 160]);
    
    
    % System parameters
    nx = PS.nx;
    nu = PS.nu;
    nw = PS.nw;
    N = PS.N;
    ScriptA = PS.ScriptA;
    ScriptB = PS.ScriptB;
    ScriptD = PS.ScriptD;
    mu0 = PS.mu0;
    muf = PS.muf;
    Sigma0 = PS.Sigma0;
    Sigmaf = PS.Sigmaf;
    A = PS.A;
    B = PS.B;
    D = PS.D;
    
    % Optimal mean motion
    EX = ScriptA * mu0 + ScriptB * V;
    nominal = value(EX);
    
    % Set random number generator
    rng(0);
    
    % Run Monte Carlo's to generate optimal trajectories
    for mc = 1:MCnum
       % Generate initial condition from normal distribution
       x0_MC = mu0 + sqrtm(Sigma0) * randn(nx,1);
       
       % Initiate variable histories
       U = zeros(nu,N); % optimal control history
       x_MC = zeros(nx,N+1); % optimal state history
       y_MC = zeros(nx,N+1); % optimal deviation history
       
       % Store initial states
       x_MC(:,1) = x0_MC;
       y_MC(:,1) = x0_MC - mu0;
       
       % Apply dynamics and store history
       for k = 1:N
           % Compute optimal control at time tk
           U(:,k) =  V((k-1)*nu+1:k*nu) + K((k-1)*nu+1:k*nu,(k-1)*nx+1:k*nx) * y_MC(:,k);
           
           % Propagate dynamics to time tk
           w = randn(nw,1);
           x_MC(:,k+1) = A * x_MC(:,k) + B * U(:,k) + D * w;
           y_MC(:,k+1) = A * y_MC(:,k) + D * w;
       end       
       
       % Plot current Monte Carlo run
       figure(fignum); hold on;
       plot(x_MC(1,:),x_MC(2,:),'color',grey);
       hold off;
    end
    
    % Set confidence on covariance ellipses
    rl = 0.9973; % 3-sigma   
    
    figure(fignum);hold on;
    
    % Plot mean motion
    plot(nominal(1:nx:end),nominal(2:nx:end),'+k','LineWidth',2,'MarkerSize',10);
    plot(nominal(1:nx:end),nominal(2:nx:end),'k','LineWidth',1);
    
    % Plot initial and terminal covariance ellipses
    h0 = error_ellipse('C',Sigma0(1:2,1:2),'mu',mu0(1:2),'conf',rl,'style','--k');
    set(h0,'linewidth',3);
    hf = error_ellipse('C',Sigmaf(1:2,1:2),'mu',muf(1:2),'conf',rl,'style','k');
    set(hf,'linewidth',2);
    
    % Useful matrices
    IplusBK = eye((N+1)*nx) + ScriptB * K;
    SigmaY = ScriptA * Sigma0 * ScriptA' + ScriptD * ScriptD';
    
    % Full state covariance matrix
    SX = IplusBK * SigmaY * IplusBK';
    % Plot history of covariance ellipses
    for k = 0:N
        Ek = [zeros(nx,k*nx) eye(nx) zeros(nx,(PS.N-k)*nx)];
        % Extract state mean and covariance at time tk
        muk = Ek * nominal;
        SXk = Ek * SX * Ek';
        
        % Plot covariance ellipse
        h = error_ellipse('C',SXk(1:2,1:2),'mu',muk(1:2),'conf',rl,'style','k');
        set(h,'linewidth',1);
    end

    xlabel('$x$ (m)','Interpreter','latex','FontSize',15);
    ylabel('$y$ (m)','Interpreter','latex','FontSize',15);
    set(gca,'FontSize',15);
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
    
    %{
    axes('position',[.65 .15 .15 .2])
    box on % put box around new pair of axes
    hold on
    rewind = 3;
    for mc = 1:MCnum
       x0_MC = PS.mu0 + sqrtm(PS.Sigma0)*randn(nx,1);
       U = zeros(nu,N);
       x_MC = zeros(nx,N+1);
       y_MC = zeros(nx,N+1);
       x_MC(:,1) = x0_MC;
       y_MC(:,1) = x0_MC - PS.mu0;
       for k = 1:N
           U(:,k) =  V((k-1)*nu+1:k*nu) + K((k-1)*nu+1:k*nu,(k-1)*nx+1:k*nx) * y_MC(:,k);
           w = randn(nw,1);
           x_MC(:,k+1) = PS.A*x_MC(:,k) + PS.B*U(:,k) + PS.D*w;
           y_MC(:,k+1) = PS.A*y_MC(:,k) + PS.D*w;
       end       
       plot(x_MC(1,end-rewind:end),x_MC(2,end-rewind:end),'color',grey);
    end
    
    nominal_x = nominal(1:nx:end);
    nominal_y = nominal(2:nx:end);
    plot(nominal_x(end-rewind-1:end),nominal_y(end-rewind-1:end),'+k','LineWidth',2,'MarkerSize',10);
    plot(nominal_x(end-rewind-1:end),nominal_y(end-rewind-1:end),'k','LineWidth',1);
    
    hf = error_ellipse('C',PS.Sigmaf(1:2,1:2),'mu',PS.muf(1:2),'conf',rl,'style','k');
    set(hf,'linewidth',2);
    
    for k = N-rewind:N
        Ek = [zeros(nx,k*nx) eye(nx) zeros(nx,(PS.N-k)*nx)];
        muk = Ek*nominal;
        SXk = Ek*SX*Ek';
        h = error_ellipse('C',SXk(1:2,1:2),'mu',muk(1:2),'conf',rl,'style','k');
        set(h,'linewidth',1);
    end
    
    %set(gca,'visible','off')
    set(gca,'XTick',[], 'YTick', [])
    %}
end