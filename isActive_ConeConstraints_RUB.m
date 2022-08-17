function [bool,diff,deltaTrue_RUB_max] = isActive_ConeConstraints_RUB(V,K,f,PS)
eps = 1E-06;
N = PS.N;

deltax = PS.deltax;

deltaTrue_RUB = computeTrueDeltaC_ConeConstraints_RUB(PS,V,K,f);  
deltaTrue_RUB_max = max(max(deltaTrue_RUB));
deltaTrue_RUB_max = deltaTrue_RUB_max(:);

diff = zeros(N,1);
bool = false(N,1);

for k = 1:N
    diff(k) = abs(deltax(k) - deltaTrue_RUB_max(k));
    if  diff(k) < eps
        bool(k) = true;
    end
end

%bool = any(bool);

end