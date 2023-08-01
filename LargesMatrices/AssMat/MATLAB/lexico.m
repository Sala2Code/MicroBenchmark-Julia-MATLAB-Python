function Num = lexico(i,j)
% Lexicographic Numbering
global Nxy;
Ny=Nxy(2);
Num=(i-1)*Ny+j;
end