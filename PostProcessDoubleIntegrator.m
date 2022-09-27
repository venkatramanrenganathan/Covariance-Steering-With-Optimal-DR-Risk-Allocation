function PostProcessDoubleIntegrator(fignum,MCnum,V,K,PS)
    %figure(fignum);clf;
    nx = PS.nx;
    nu = PS.nu;
    nw = PS.nw;
    N = PS.N;
    A = PS.A;
    B = PS.B;
    D = PS.D;
    mu0 = PS.mu0;
    Sigma0 = PS.Sigma0;
    muf = PS.muf;
    Sigmaf = PS.Sigmaf;
    
    grey = [0.8 0.8 0.8];
    darkgrey = [0.6 0.6 0.6];
    orange = [0.850 0.325 0.098];
    
    EX = PS.ScriptA*PS.mu0 + PS.ScriptB*V;
    
    rng(0);
    for mc = 1:MCnum
       x0_MC = mvlaprnd(nx,PS.mu0,PS.Sigma0,1); 
       %x0_MC = PS.mu0 + sqrtm(PS.Sigma0)*randn(nx,1);
       U = zeros(nu,N);
       x_MC = zeros(nx,N+1);
       y_MC = zeros(nx,N+1);
       x_MC(:,1) = x0_MC;
       y_MC(:,1) = x0_MC - mu0;
       % Create disturbance history
       wMean = zeros(nw, 1);
       wCov  = eye(nw);
       wData = mvlaprnd(nw,wMean,wCov,N);
       for k = 1:N
           U(:,k) =  V((k-1)*nu+1:k*nu) + K((k-1)*nu+1:k*nu,(k-1)*nx+1:k*nx) * y_MC(:,k);
           w = wData(:,k); % randn(nw,1);
           x_MC(:,k+1) = A * x_MC(:,k) + B * U(:,k) + D * w;
           y_MC(:,k+1) = A * y_MC(:,k) + D * w;
       end       
       
       figure(fignum); hold on;
       plot(x_MC(1,:),x_MC(2,:),'color',orange); 
       a = findobj(gcf, 'type', 'axes');
       h = findobj(gcf, 'type', 'line');
       set(h, 'linewidth', 4);
       set(a, 'linewidth', 4);
       set(a, 'FontSize', 40);
       hold off;
    end
    nominal = value(EX);
    rl = 0.9973;% 3 Sigma   
    
    figure(fignum);
    hold on;
    plot(nominal(1:nx:end),nominal(2:nx:end),'+k','LineWidth',2,'MarkerSize',10);
    plot(nominal(1:nx:end),nominal(2:nx:end),'k','LineWidth',1);
    h0 = error_ellipse('C',Sigma0(1:2,1:2),'mu',mu0(1:2),'conf',rl,'style','--k');
    set(h0,'LineWidth',3);
    hf = error_ellipse('C',Sigmaf(1:2,1:2),'mu',muf(1:2),'conf',rl,'style','k');
    set(hf,'LineWidth',2);
    IplusBK = eye((N+1)*nx)+PS.ScriptB*K;
    SigmaY = PS.ScriptA*PS.Sigma0*PS.ScriptA'+PS.ScriptD*PS.ScriptD';
    SX = IplusBK*SigmaY*IplusBK';
    for k = 0:N
        Ek = [zeros(nx,k*nx) eye(nx) zeros(nx,(N-k)*nx)];
        muk = Ek * nominal;
        SXk = Ek * SX * Ek';
        if k == N
            h = error_ellipse('C',SXk(1:2,1:2),'mu',muk(1:2),'conf',rl,'style','k');
            set(h,'LineWidth',1);
            h = error_ellipse('C',SXk(1:2,1:2),'mu',muk(1:2),'conf',rl,'style','--k');
            set(h,'LineWidth',3);
        else
            h = error_ellipse('C',SXk(1:2,1:2),'mu',muk(1:2),'conf',rl,'style','k');
            set(h,'LineWidth',1);
        end
    end
    
    xl = linspace(-20,0.4);
    %plot(xl, 0.2*(-xl + PS.stateCC_offset), '-b','LineWidth',1);
    %plot(xl, 0.2*(xl + PS.stateCC_offset), '-b','LineWidth',1);
    plot(xl, 0.2*(-xl + 1), '-b','LineWidth',1);
    plot(xl, 0.2*(xl - 1), '-b','LineWidth',1);
    %ineqplot('y > 0.2*x-0.2',[-20 10 -10 10],darkgrey);
    %ineqplot('y < -0.2*x+0.2',[-20 10 -10 10],darkgrey);
    xlabel('$x$ (m)','Interpreter','latex','FontSize',15);
    ylabel('$y$ (m)','Interpreter','latex','FontSize',15);
    grid on;
    a = findobj(gcf, 'type', 'axes');
    h = findobj(gcf, 'type', 'line');
    set(h, 'linewidth', 4);
    set(a, 'linewidth', 4);
    set(a, 'FontSize', 40);
    xlim([-20, 5]);
    ylim([-4, 4]);
    hold off;
end