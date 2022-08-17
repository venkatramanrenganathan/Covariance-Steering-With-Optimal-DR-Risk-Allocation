function PlotResultsRPO(Vstar,Kstar,PS)
    figure(2); 
    subplot(1,2,1); hold on; grid on;
    plot(PS.dt * [1:PS.N],PS.deltax(2,:),'.k','MarkerSize',20);
    plot(PS.dt * [1:PS.N],PS.deltax(3,:),'ob','MarkerSize',10,'LineWidth',1);
    legend({'$\delta_r$','$\delta_u$'},'Interpreter','latex','FontSize',20);
    xlabel('Time (s)','Interpreter','latex','FontSize',15);
    ylabel('Allocated Risk','FontSize',15);
    set(gca,'FontSize',15);
    subplot(1,2,2); hold on; grid on;
    plot(PS.dt * [1:PS.N],PS.deltabar(2,:),'.k','MarkerSize',20);
    plot(PS.dt * [1:PS.N],PS.deltabar(3,:),'ob','MarkerSize',10,'LineWidth',1);
    legend({'$\bar{\delta}_r$','$\bar{\delta}_u$'},'Interpreter','latex','FontSize',20);
    xlabel('Time (s)','Interpreter','latex','FontSize',15);
    ylabel('True Risk','FontSize',15);
    set(gca,'FontSize',15);

    PostProcess(3,100,Vstar,Kstar,PS);
    PostProcess3D(4,100,Vstar,Kstar,PS,1);
    PostProcess3D(7,100,Vstar,Kstar,PS,0);
    xlim([-0.2 0.2]); ylim([-0.2 0.1]); zlim([-0.2 0.2]);
end