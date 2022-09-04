function PostProcess3D(fignum,MCnum,V,K,PS,bool)
    warning('off','all')
    %% Plot Polyehdral State Space
    bound = 150;
    x = linspace(-bound,0,100000);
    y = x;
    z = x;
    figure(fignum); 
    plot3(x,y+PS.stateCC_offset,z,'k','LineWidth',2); hold on; grid on;
    plot3(-x,y+PS.stateCC_offset,z,'k','LineWidth',2);
    plot3(x,y+PS.stateCC_offset,-z,'k','LineWidth',2);
    plot3(-x,y+PS.stateCC_offset,-z,'k','LineWidth',2);
    plot3(bound*ones(1,100),-(bound-PS.stateCC_offset)*ones(1,100),linspace(-bound,bound),'k','LineWidth',2);
    plot3(-bound*ones(1,100),-(bound-PS.stateCC_offset)*ones(1,100),linspace(-bound,bound),'k','LineWidth',2);
    plot3(linspace(-bound,bound),-(bound-PS.stateCC_offset)*ones(1,100),-bound*ones(1,100),'k','LineWidth',2);
    plot3(linspace(-bound,bound),-(bound-PS.stateCC_offset)*ones(1,100),bound*ones(1,100),'k','LineWidth',2);
%     plot3(x,y+.05,zeros(1,100),'--b','LineWidth',0.5);
%     plot3(-x,y+.05,zeros(1,100),'--b','LineWidth',0.5);
%     plot3(linspace(-1.25,1.25),-1.2*ones(1,100),zeros(1,100),'--b','LineWidth',0.5);
    v1 = [0 0];
    v2 = [-bound -(bound-PS.stateCC_offset)];
    v3 = [bound -(bound-PS.stateCC_offset)];
    a1 = [v1(1) v2(1) v3(1)];
    a2 = [v1(2) v2(2) v3(2)];
    patch(a1,a2,'b');
    plot3(0,0,0,'.r','MarkerSize',20);
    xlabel('$x$ (m)','Interpreter','latex','FontSize',15); 
    ylabel('$y$ (m)','Interpreter','latex','FontSize',15); 
    zlabel('$z$ (m)','Interpreter','latex','FontSize',15);
    alpha(0.3);

    %% Plot optimal trajectories and covariance ellipsoids
    grey    = [0.7 0.7 0.7];
    nx      = PS.nx;
    nu      = PS.nu;
    nw      = PS.nw;
    N       = PS.N;
    
    mu0     = PS.mu0;
    muf     = PS.muf;
    Sigma0  = PS.Sigma0;
    Sigmaf  = PS.Sigmaf;
    
    ScriptA = PS.ScriptA;
    ScriptB = PS.ScriptB;
    ScriptD = PS.ScriptD;
    
    A       = PS.A;
    B       = PS.B;
    D       = PS.D;
    
    rng(0);
    
    U = zeros(nu, N, MCnum);
    xhist = zeros(nx, N + 1, MCnum);
    
    for mc = 1 : MCnum
        
       x0_MC     = mu0 + sqrtm(Sigma0) * randn(nx,1);
       x_MC      = zeros(nx, N+1);
       y_MC      = zeros(nx, N+1);
       x_MC(:,1) = x0_MC;
       y_MC(:,1) = x0_MC - mu0;
       wMean = zeros(nw, 1);
       wCov  = eye(nw);
       wData = mvlaprnd(nw,wMean,wCov,N);
       
       for k = 1 : N
           
           U(:, k, mc) = V((k-1)*nu+1:k*nu) + K((k-1)*nu+1:k*nu,(k-1)*nx+1:k*nx) * y_MC(:,k);
           wk          = wData(:, k); % randn(nw,1);
           x_MC(:,k+1) = A * x_MC(:,k) + B * U(:, k, mc) + D * wk;
           y_MC(:,k+1) = A * y_MC(:,k) + D * wk;
       
       end       
       
       xhist(:, :, mc) = x_MC;
       plot3(x_MC(1, :), x_MC(2, :), x_MC(3, :), 'color', grey);
    end
    
    EX         = ScriptA * mu0 + ScriptB * V;
    nominal    = value(EX);
    px_nominal = nominal(1 : nx : end);
    py_nominal = nominal(2 : nx : end);
    pz_nominal = nominal(3 : nx : end);
    rl         = 0.9973; % 3 Sigma   

    plot3(px_nominal, py_nominal, pz_nominal, '+k', 'LineWidth', 2, 'MarkerSize',10);
    plot3(px_nominal, py_nominal, pz_nominal , 'k', 'LineWidth', 1);
    h0 = error_ellipse('C', Sigma0(1 : 3, 1 : 3), 'mu', mu0(1 : 3), 'conf', rl, 'style', 'r');
    shading interp;
    %set(h0,'linewidth',0.5);
    hf = error_ellipse('C',PS.Sigmaf(1:3,1:3),'mu',PS.muf(1:3),'conf',rl,'style','r');
    shading interp;
    %set(hf,'linewidth',0.5);
    
    IplusBK = eye((N + 1) * nx) + ScriptB * K;
    SigmaY  = ScriptA * Sigma0 * ScriptA' + ScriptD * ScriptD';
    SX      = IplusBK * SigmaY * IplusBK';
    stdhist = zeros(nx, N+1);
    for k = 0 : N
        
        Ek  = [zeros(nx, k * nx) eye(nx) zeros(nx, (N - k) * nx)];
        muk = Ek * nominal;
        SXk = Ek * SX * Ek';
        stdhist(:,k+1) = sqrt(diag(SXk));
        h = error_ellipse('C', SXk(1 : 3, 1 : 3), 'mu', muk(1 : 3), 'conf', rl, 'style', 'k');
        shading interp;
        %set(h,'linewidth',0.5);
    end
    
    if bool
        %% Plot State Histories and Standard Deviations
        grey = [0.8 0.8 0.8];
        figure(fignum + 1);
        timeStepVec_x = linspace(0, N, N + 1);
        timeVec_x = PS.dt .* timeStepVec_x;
        
        subplot(3,2,1); hold on; grid on;
        for mc = 1:MCnum
            plot(timeVec_x,xhist(1,:,mc),'color',grey);
        end
        plot(timeVec_x,nominal(1:nx:end),'+k','LineWidth',2,'MarkerSize',10);
        plot(timeVec_x,nominal(1:nx:end),'k','LineWidth',1);
        shadedErrorBar(timeVec_x,nominal(1:nx:end),3*stdhist(1,:),'lineprops','-k','transparent',1,'patchSaturation',0.1);
        xlabel('Time (s)','FontSize',13);
        ylabel('$x$ (m)','Interpreter','latex','FontSize',13);
        set(gca,'FontSize',11);
        
%         subplot(3,2,2); hold on; grid on;
%         plot(timeVec_x,stdhist(1,:),'-k','LineWidth',2);
%         xlabel('Time (s)','FontSize',13);
%         ylabel('$\sigma_x$ (m)', 'Interpreter','latex','FontSize',13);
%         set(gca,'FontSize',11);

        subplot(3,2,3); hold on; grid on;
        for mc = 1:MCnum
            plot(timeVec_x,xhist(2,:,mc),'color',grey);
        end
        plot(timeVec_x,nominal(2:nx:end),'+k','LineWidth',2,'MarkerSize',10);
        plot(timeVec_x,nominal(2:nx:end),'k','LineWidth',1);
        shadedErrorBar(timeVec_x,nominal(2:nx:end),3*stdhist(2,:),'lineprops','-k','transparent',1,'patchSaturation',0.1);
        xlabel('Time (s)','FontSize',13);
        ylabel('$y$ (m)','Interpreter','latex','FontSize',13);
        set(gca,'FontSize',11);
        
%         subplot(3,2,4); hold on; grid on;
%         plot(timeVec_x,stdhist(2,:),'-k','LineWidth',2);
%         xlabel('Time (s)','FontSize',13);
%         ylabel('$\sigma_y$ (m)', 'Interpreter','latex','FontSize',13);
%         set(gca,'FontSize',11);

        subplot(3,2,5); hold on; grid on;
        for mc = 1:MCnum
            plot(timeVec_x,xhist(3,:,mc),'color',grey);
        end
        plot(timeVec_x,nominal(3:nx:end),'+k','LineWidth',2,'MarkerSize',10);
        plot(timeVec_x,nominal(3:nx:end),'k','LineWidth',1);
        shadedErrorBar(timeVec_x,nominal(3:nx:end),3*stdhist(3,:),'lineprops','-k','transparent',1,'patchSaturation',0.1);
        xlabel('Time (s)','FontSize',13);
        ylabel('$z$ (m)','Interpreter','latex','FontSize',13);
        set(gca,'FontSize',11);
        
%         subplot(3,2,6); hold on; grid on;
%         plot(timeVec_x,stdhist(3,:),'-k','LineWidth',2);
%         xlabel('Time (s)','FontSize',13);
%         ylabel('$\sigma_z$ (m)', 'Interpreter','latex','FontSize',13);
%         set(gca,'FontSize',11);

        %% Plot optimal controls
        timeStepVec_u = 0:N-1;
        timeVec_u = PS.dt .* timeStepVec_u;
        
        SigmaUU = K * SigmaY * K';
        SigmaUU_hist = zeros(nu,nu,N);
        stdUU_hist = zeros(nu,N);
        for k = 0:N-1
            Fk = [zeros(nu,k*nu) eye(nu) zeros(nu,(N-k-1)*nu)];
            SigmaUU_hist(:,:,k+1) = Fk * SigmaUU * Fk';
            for i = 1:nu
                stdUU_hist(i,k+1) = sqrt(SigmaUU_hist(i,i,k+1));
            end
        end
    
    
        figure(fignum+1); 
        for k = 1:nu
            subplot(3,2,2*k); hold on; grid on;
            xlabel('Time (s)','FontSize',13);
            ylabel(sprintf('Thruster %i (N)',k),'Interpreter','latex','FontSize',13);
            for mc = 1:MCnum
                plot(timeVec_u,U(k,:,mc),'color',grey);
            end
            plot(timeVec_u,V(k:nu:end),'+k','LineWidth',2,'MarkerSize',10);
            plot(timeVec_u,V(k:nu:end),'k','LineWidth',1);
            shadedErrorBar(timeVec_u,V(k:nu:end),3*stdUU_hist(k,:),'lineprops','-k','transparent',1,'patchSaturation',0.1);
            plot(timeVec_u,ones(1,N)*PS.umax,'--k','LineWidth',1);
            plot(timeVec_u,-ones(1,N)*PS.umax,'--k','LineWidth',1);
            set(gca,'FontSize',11);
            ylim([-(PS.umax+PS.umax/4),PS.umax+PS.umax/4]);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

    end
end