clear all;

global beta rTol aTol;

beta = 3;
aTol = 10^(-10);
rTol = 10^(-10);
NpointsX = 15;
NpointsY = 15;
    
X = linspace(0, 1, NpointsX); 
Y = linspace(0, 1,  NpointsY);

solTheorique = zeros(NpointsY, NpointsX);
solPropagee = zeros(NpointsY, NpointsX);
solError  = zeros(NpointsY, NpointsX);

% partie peu chère : valeurs théoriques 
for i = NpointsX-1:-1:2
    for j = NpointsY-1:-1:2
         xIci = X(i);
         yIci = Y(j);
         % changement i,j devient j,i
        solTheorique(j,i) = Sol(xIci, yIci);
        [timeFoot, foot] = propagation([xIci, yIci]);
        solPropagee(j,i) = Sol(foot(1), foot(2));
        solError(j,i) = abs(solTheorique(j,i) - solPropagee(j,i));
    end
end

