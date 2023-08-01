function MatSol=Sol(x,y)
global beta
t2 = x ^ 2;
t6 = cos(pi * y);
MatSol = sin(pi * x + t6 * (t2 - x) * beta);
