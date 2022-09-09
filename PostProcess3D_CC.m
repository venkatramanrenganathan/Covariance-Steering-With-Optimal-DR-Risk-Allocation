 function PostProcess3D_CC(fignum,MCnum,V,K,PS,bool)
    warning('off','all')
    %% Plot Conical State Space
    figure(fignum); 
    
    h = 200;
    r = h * tand(PS.coneAngle);
    Np = 1000;
    
    u = linspace(0,h,Np) ;
    th = linspace(0,2*pi,Np) ;
    
    [U,Th] = meshgrid(u,th);
    X = ((h - U) / h) .* r .* cos(Th);
    Y = -U + h - PS.d;
    Z = ((h - U) / h) .* r .* sin(Th);
    
    Rx = [cosd(PS.psi) sind(PS.psi);
          -sind(PS.psi) cosd(PS.psi)];
    for ii = 1 : Np
        for jj = 1 : Np
            yz_new = Rx * [Y(ii,jj); Z(ii,jj)];
            Y(ii,jj) = yz_new(1);
            Z(ii,jj) = yz_new(2);
        end
    end
            
    hSurface = surf(X,Y,Z); 
    set(hSurface,'FaceColor',[0 0 1], ...
      'FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor','none');
    hold on;
  
    %% Plot optimal trajectories and covariance ellipsoids
    % Color for plotting
    grey = [0.7 0.7 0.7];
    darkgrey = [0.5 0.5 0.5];
    
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
    rng(2);
    
    % Initial optimal state and control histories
    U = zeros(nu,N,MCnum);
    xhist = zeros(nx,N+1,MCnum);
    for mc = 1:MCnum
       x0_MC = mu0 + sqrtm(Sigma0) * randn(nx,1);
       x_MC = zeros(nx,N+1);
       y_MC = zeros(nx,N+1);
       x_MC(:,1) = x0_MC;
       y_MC(:,1) = x0_MC - mu0;
       wMean = zeros(nw, 1);
       wCov  = eye(nw);
       wData = mvlaprnd(nw,wMean,wCov,N);
       for k = 1:N
           U(:,k,mc) =  V((k-1)*nu+1:k*nu) + K((k-1)*nu+1:k*nu,(k-1)*nx+1:k*nx) * y_MC(:,k);
           w = wData(:, k); % randn(nw,1);
           x_MC(:,k+1) = A * x_MC(:,k) + B * U(:,k,mc) + D * w;
           y_MC(:,k+1) = A * y_MC(:,k) + D * w;
       end       
       xhist(:,:,mc) = x_MC;
       plot3(x_MC(1,:),x_MC(2,:),x_MC(3,:),'color',darkgrey);
    end
    
    % Set confidence on covariance ellipsoids
    rl = 0.9973; % 3-sigma   
    
    % Plot optimal mean trajectory
    plot3(nominal(1:nx:end),nominal(2:nx:end),nominal(3:nx:end),'+k','LineWidth',2,'MarkerSize',10);
    plot3(nominal(1:nx:end),nominal(2:nx:end),nominal(3:nx:end),'k','LineWidth',1);
    
    % Plot initial and terminal ellipsoids
    h0 = error_ellipse('C',PS.Sigma0(1:3,1:3),'mu',PS.mu0(1:3),'conf',rl,'style','r');
    shading interp;
    %set(h0,'linewidth',0.5);
    hf = error_ellipse('C',PS.Sigmaf(1:3,1:3),'mu',PS.muf(1:3),'conf',rl,'style','r');
    shading interp;
    %set(hf,'linewidth',0.5);
    
    % Useful matrices
    IplusBK = eye((N+1)*nx) + ScriptB * K;
    SigmaY = ScriptA * Sigma0 * ScriptA' + ScriptD * ScriptD';
    
    % Full state covariance
    SX = IplusBK * SigmaY * IplusBK';
    
    % Initialize state standard deviation histories
    stdhist = zeros(nx,N+1);
    for k = 0:N
        Ek = [zeros(nx,k*nx) eye(nx) zeros(nx,(PS.N-k)*nx)];
        
        % Compute state mean and covariance at time tk
        muk = Ek * nominal;
        SXk = Ek * SX * Ek';
        
        % Store state standard deviation at time tk
        stdhist(:,k+1) = sqrt(diag(SXk));
        
        % Plot state covariance ellipsoids at time tk
        h = error_ellipse('C',SXk(1:3,1:3),'mu',muk(1:3),'conf',rl,'style','k');
        shading interp;
        %set(h,'linewidth',0.5);
    end
    
    if bool
        %% Plot State Histories and Standard Deviations
        grey = [0.8 0.8 0.8];
        figure(fignum+1);
        timeStepVec_x = linspace(0,N,N+1);
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
            %ylim([-(PS.umax+PS.umax/4),PS.umax+PS.umax/4]);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

    end
end