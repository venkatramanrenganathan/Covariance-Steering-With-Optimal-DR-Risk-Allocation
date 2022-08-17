function PostProcess(fignum,MCnum,V,K,PS)
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
    
    EX = PS.ScriptA*PS.mu0 + PS.ScriptB*V;
    
    rng(0);
    for mc = 1:MCnum
       x0_MC = PS.mu0 + sqrtm(PS.Sigma0)*randn(nx,1);
       U = zeros(nu,N);
       x_MC = zeros(nx,N+1);
       y_MC = zeros(nx,N+1);
       x_MC(:,1) = x0_MC;
       y_MC(:,1) = x0_MC - mu0;
       for k = 1:N
           U(:,k) =  V((k-1)*nu+1:k*nu) + K((k-1)*nu+1:k*nu,(k-1)*nx+1:k*nx) * y_MC(:,k);
           w = randn(nw,1);
           x_MC(:,k+1) = A * x_MC(:,k) + B * U(:,k) + D * w;
           y_MC(:,k+1) = A * y_MC(:,k) + D * w;
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
%             h = error_ellipse('C',SXk(1:2,1:2),'mu',muk(1:2),'conf',rl,'style','k');
%             set(h,'LineWidth',1);
    end
    
    xl = linspace(-20,0);
    xr = linspace(0,140);
    plot(xr, -xr + PS.stateCC_offset, '--k','LineWidth',1);
    plot(xl, xl + PS.stateCC_offset, '--k','LineWidth',1);
    ineqplot('y > -x+5',[-20 140 -140 10],darkgrey);
    ineqplot('y > x+5',[-20 140 -140 10],darkgrey);
    grid on;
    xlabel('$x$ (m)','Interpreter','latex','FontSize',15);
    ylabel('$y$ (m)','Interpreter','latex','FontSize',15);
    hold off;
end