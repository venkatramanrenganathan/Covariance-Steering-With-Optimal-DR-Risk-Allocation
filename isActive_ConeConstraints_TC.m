function [bool,diff,deltaTrue_TC_max] = isActive_ConeConstraints_TC(V,f,t,PS)
eps = 1E-06;
N = PS.N;

deltax = PS.deltax;

deltaTrue_TC = computeTrueDeltaC_ConeConstraints_ThreeCut(PS,V,f,t); 
deltaTrue_TC_max = max(max(deltaTrue_TC));
deltaTrue_TC_max = deltaTrue_TC_max(:);

diff = zeros(N,1);
bool = false(N,1);

for k = 1:N
    diff(k) = abs(deltax(k) - deltaTrue_TC_max(k));
    if  diff(k) < eps
        bool(k) = true;
    end
end

%bool = any(bool);

end