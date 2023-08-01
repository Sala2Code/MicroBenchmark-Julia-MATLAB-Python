for j=1:Ny
   Numl = lexico(1,j);
   Numc = lexico(2,j);
   Mat(Numl,:) = 0;
   Mat(Numl,Numl) = 1/(Dx*Dx);
   %Mat(Numl,Numc) = 1/(Dx*Dx);
   RhsV(Numl) = SolM(1,j)*Mat(Numl,Numl);
   
   
   Numl = lexico(Nx,j);
   Numc = lexico(Nx-1,j);
   Mat(Numl,:) = 0;
   Mat(Numl,Numl) = 1/(Dx*Dx);
   %Mat(Numl,Numc) = 1/(Dx*Dx);
   RhsV(Numl) = SolM(Nx,j)*Mat(Numl,Numl);
end
% 
for i=2:Nx-1
   xi = x(i);
   Numl = lexico(i,1);
   yj = (y(1)+y(2))*0.5;
   RhsV(Numl) = yTraceBeta(xi,yj,Beta,epsilon,nPeriod)/Dy; % negative outward normal
   Numl = lexico(i,Ny);
   yj = (y(Ny)+y(Ny-1))*0.5;
   RhsV(Numl) = -yTraceBeta(xi,yj,Beta,epsilon,nPeriod)/Dy;
end

% Dirichlet Boundary 
% for i=2:Nx-1
%    Numl = lexico(i,1);
%    Mat(Numl,:) = 0; Mat(Numl,Numl) = 1/(Dy*Dy);
%    RhsV(Numl) = SolM(i,1)*Mat(Numl,Numl);
%    Numl = lexico(i,Ny);
%    Mat(Numl,:) = 0; Mat(Numl,Numl) = 1/(Dy*Dy);
%    RhsV(Numl) = SolM(i,Ny)*Mat(Numl,Numl);
% end