clear all
%% Parameters
global epsilon;
epsilon = 1.e-2;
global nPeriod;
nPeriod = 1;
global Beta;
Beta = 1;
k = 50; % nodes
%% Grid Definition
global Dx Dy

Ly = 1; 
Lx = Ly; 
Nx = k;
Ny = k;

Dx=Lx/(Nx-2); x=(-Dx/2:Dx:Lx+Dx/2);
Dy=Ly/(Ny-2); y=(-Dy/2:Dy:Ly+Dy/2); % Match the physical boundary at mid points

global Nxy
Nxy = [Nx; Ny];
NL = Nx*Ny; % Number of lines of the linear system

t1 = cputime;
%% Problem setup
X = zeros(Nx,Ny); Y = zeros(Nx,Ny);
Matb1=zeros(Nx,Ny); Matb2=zeros(Nx,Ny);
SolM = zeros(Nx,Ny);
dxSolM = zeros(Nx,Ny); dySolM = dxSolM;
dxxSolM = dxSolM; dyySolM = dySolM; dxySolM = dxSolM;
RhsM = zeros(Nx,Ny);
bGradSolM = zeros(Nx,Ny);
SolV = zeros(NL,1) ; RhsV = zeros(NL,1);
SolV1= zeros(NL,1);

B1xdy = zeros(Nx,Ny); B2xdy = zeros(Nx,Ny);
B1xyd = zeros(Nx,Ny); B2xyd = zeros(Nx,Ny);

dxBx = zeros(Nx,Ny); dyBx = zeros(Nx,Ny);
dxBy = zeros(Nx,Ny); dyBy = zeros(Nx,Ny);

dxdSolM = zeros(Nx,Ny); dydSolM = zeros(Nx,Ny);
BxdxU1d = zeros(Nx,Ny); BydyU1d = zeros(Nx,Ny);
SolM1 = zeros(Nx,Ny); 
B1xdyd =zeros(Nx-1,Ny-1); B2xdyd=B1xdyd;
for i=1:Nx-1
    xd = (x(i)+x(i+1))/2;
    for j=1:Ny-1
        yd = ( y(j)+y(j+1))/2;
            
        B1xdy(i,j) = BxBeta(xd,y(j),Beta);
        B2xdy(i,j) = ByBeta(xd,y(j),Beta);
        B1xyd(i,j) = BxBeta(x(i),yd,Beta);
        B2xyd(i,j) = ByBeta(x(i),yd,Beta);
        B1xdyd(i,j) = BxBeta(xd,yd,Beta);
        B2xdyd(i,j) = ByBeta(xd,yd,Beta);
    end
end


for i=1:Nx
    for j=1:Ny
        Numl=lexico(i,j);
        X(i,j) = x(i);        
        Y(i,j) = y(j);
        % Magnetic field
        Matb1(i,j)=BxBeta(x(i),y(j),Beta);
        Matb2(i,j)=ByBeta(x(i),y(j),Beta);
                
        SolM(i,j)=SolBeta(x(i),y(j),Beta,epsilon,nPeriod);
        LapParaSolM(i,j) = LapParaSolBeta(x(i),y(j),Beta,epsilon,nPeriod);
        LapOrthoSolM(i,j) = LapOrthoSolBeta(x(i),y(j),Beta,epsilon,nPeriod);
        SolV(Numl)=SolM(i,j);
        RhsM(i,j)=RhsBeta(x(i),y(j),Beta,epsilon,nPeriod)/epsilon;
        RhsV(Numl)=RhsM(i,j);
        bGradSolM(i,j)=bGradSolBeta(x(i),y(j),Beta,epsilon,nPeriod);
        
    end
end


RhsVbis = RhsV;
Mat11 = Matb1.*Matb1;
Mat12 = Matb1.*Matb2;
Mat22 = Matb2.*Matb2;
Mat21 = Mat12;

Matb11xdy = B1xdy.* B1xdy; Matb12xdy = B1xdy.* B2xdy;
Matb22xyd = B2xyd.* B2xyd; Matb21xyd = B1xyd.* B2xyd;
Mat21xdyd = B1xdyd .* B2xdyd;
MatOnes = ones(Nx,Ny); MatZeros = zeros(Nx,Ny);

% Flux Sym"
LapPara = LapParaSym(Mat11,Mat12,Mat21,Mat22);
Lap    = LapParaSym(MatOnes,MatZeros,MatZeros,MatOnes);

LapOrtho = Lap - LapPara;
Mat = LapOrtho + LapPara/epsilon;
SetBoundariesBeta;

%% Solve System 
 Approx = Mat\RhsV;
 ApproxV = reshape(Approx,Ny,Nx)';
 Error = ApproxV - SolM;

% disp(max(abs(Error(:))));
