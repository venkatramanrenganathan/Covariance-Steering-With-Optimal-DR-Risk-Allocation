function PostProcessCC(fignum,MCnum,V,K,PS)
    %figure(fignum);clf;
    grey = [0.7 0.7 0.7];
    nx = PS.nx;
    nu = PS.nu;
    nw = PS.nw;
    N = PS.N;
    EX = PS.ScriptA*PS.mu0 + PS.ScriptB*V;
    rng(0);
    for mc = 1:MCnum
       x0_MC = mvlaprnd(nx,PS.mu0,PS.Sigma0,1);  
       % x0_MC = PS.mu0 + sqrtm(PS.Sigma0)*randn(nx,1);
       U = zeros(nu,N);
       x_MC = zeros(nx,N+1);
       y_MC = zeros(nx,N+1);
       x_MC(:,1) = x0_MC;
       y_MC(:,1) = x0_MC - PS.mu0;
       
       % Create disturbance history
       wMean = zeros(nw, 1);
       wCov  = eye(nw);
       wData = mvlaprnd(nw,wMean,wCov,N);
       
       for k = 1:N
           U(:,k) =  V((k-1)*nu+1:k*nu) + K((k-1)*nu+1:k*nu,(k-1)*nx+1:k*nx) * y_MC(:,k);
           w = wData(:, k); % randn(nw,1);
           x_MC(:,k+1) = PS.A*x_MC(:,k) + PS.B*U(:,k) + PS.D*w;
           y_MC(:,k+1) = PS.A*y_MC(:,k) + PS.D*w;
       end       
       figure(fignum); hold on;
       plot(x_MC(1,:),x_MC(2,:),'color',grey);
       hold off;
    end
    nominal = value(EX);
    rl = 0.9973;% 3 Sigma   
    
    figure(fignum);
    hold on;
    plot(nominal(1:nx:end),nominal(2:nx:end),'+k','LineWidth',2,'MarkerSize',10);
    plot(nominal(1:nx:end),nominal(2:nx:end),'k','LineWidth',1);
    h0 = error_ellipse('C',PS.Sigma0(1:2,1:2),'mu',PS.mu0(1:2),'conf',rl,'style','--k');
    set(h0,'linewidth',3);
    hf = error_ellipse('C',PS.Sigmaf(1:2,1:2),'mu',PS.muf(1:2),'conf',rl,'style','k');
    set(hf,'linewidth',2);
    
    IplusBK = eye((N+1)*nx)+PS.ScriptB*K;
    SigmaY = PS.ScriptA*PS.Sigma0*PS.ScriptA'+PS.ScriptD*PS.ScriptD';
    SX = IplusBK*SigmaY*IplusBK';
    for k = 0:N
        Ek = [zeros(nx,k*nx) eye(nx) zeros(nx,(PS.N-k)*nx)];
        muk = Ek*nominal;
        SXk = Ek*SX*Ek';
        h = error_ellipse('C',SXk(1:2,1:2),'mu',muk(1:2),'conf',rl,'style','k');
        set(h,'linewidth',1);
    end

    xlabel('$x$ (m)','Interpreter','latex','FontSize',15);
    ylabel('$y$ (m)','Interpreter','latex','FontSize',15);
    set(gca,'FontSize',15);
    
    %%%%%%%%%%%
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
end