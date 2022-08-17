figure; grid on; hold on;
JstarVec_GEO = load('JstarVec_GEO.mat');
JstarVec_TS = load('JstarVec_TS.mat');
JstarVec_poly = load('JstarVec_poly.mat');

plot(JstarVec_poly.JstarVec,'.--k','MarkerSize',15);
plot(JstarVec_GEO.JstarVec,'.--b','MarkerSize',15);
plot(JstarVec_TS.JstarVec,'.--r','MarkerSize',15);

xlabel('Iteration Number','FontSize',15);
ylabel('Optimal Cost','FontSize',15);

legend({'Polyhedron','Geometric','Two-Sided'},'FontSize',15);
set(gca,'FontSize',15);