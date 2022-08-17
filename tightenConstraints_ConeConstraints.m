function PS = tightenConstraints_ConeConstraints(PS,alpha,bool,deltaTrue_max)
N = PS.N;
deltax = PS.deltax;

for k = 1:N
    if ~bool(k)
        deltaTruek = deltaTrue_max(k);
        deltax(k) = alpha * deltax(k) + (1 - alpha) * deltaTruek;
    end
end

PS.deltax = deltax;
end