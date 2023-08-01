function [temps, pied] = propagation(Ydepart)

global rTol aTol
options = odeset('Events', @evenementRenverse,'RelTol', rTol, 'AbsTol', aTol);
intervalleT = [0 15];
[testT, testY, testTE, testYE, testIE] = ode45(@champMagNormaliseRenverse, intervalleT, Ydepart, options);

pied = testYE';
temps = testTE;
return
end
