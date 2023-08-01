using BenchmarkTools
using LinearAlgebra
using SparseArrays


include("./Calculs.jl") # BxBeta - ByBeta - SolBeta - LapPara - LapOrtho - RhsBeta - bGradSolBeta - yTraceBeta
include("./Lap.jl") #  LapParaSym - LapParasym - StackBuffer 
include("./lexico.jl") # lexico

# @inbounds sur les boucles for contenant des accès à des tableaux
# Dans Lap.jl la fonction LapParaSym retourne une matrice creuse à l'aide de views, cela diminue l'allocation, + grand est le tableau + c'est efficace

function main()
    epsilon = 1e-2
    nPeriod = 1
    beta = 1
    k = 50

    Lx = 1
    Ly = 1
    Nx, Ny = k, k

    Dx = Lx/(Nx-2)
    x= -Dx/2 : Dx : Lx+Dx/2 
    Dy = Ly/(Ny-2)
    y= -Dy/2 : Dy : Ly+Dy/2
    Nxy = @SVector [Nx,Ny]
    NL = Nx*Ny

    # Pré-allocation
    X = zeros(Nx,Ny)
    Y = zeros(Nx,Ny)
    Matb1 = zeros(Nx,Ny)
    Matb2 = zeros(Nx,Ny)
    SolM = zeros(Nx,Ny)

    dxSolM = zeros(Nx,Ny)
    dySolM = zeros(Nx,Ny)
    dxxSolM = zeros(Nx,Ny)
    dyySolM = zeros(Nx,Ny)
    dxySolM = zeros(Nx,Ny)

    RhsM = zeros(Nx,Ny)
    bGradSolM = zeros(Nx,Ny)
    SolV = zeros(NL,1)
    RhsV = zeros(NL,1)
    SolV1 = zeros(NL,1)

    B1xdy = zeros(Nx,Ny)
    B2xdy = zeros(Nx,Ny)
    B1xyd = zeros(Nx,Ny)
    B2xyd = zeros(Nx,Ny)

    dxBx = zeros(Nx,Ny)
    dyBx = zeros(Nx,Ny)

    dxBy = zeros(Nx,Ny)
    dyBy = zeros(Nx,Ny)

    dxdSolM = zeros(Nx,Ny)
    dydSolM = zeros(Nx,Ny)
    BxdxU1d = zeros(Nx,Ny)
    BydyU1d = zeros(Nx,Ny)

    SolM1 = zeros(Nx,Ny)
    B1xdyd = zeros(Nx-1, Ny-1)
    B2xdyd = zeros(Nx-1, Ny-1)

    LapParaSolM = zeros(Nx,Ny)
    LapOrthoSolM = zeros(Nx,Ny)

    @inbounds  for i in 1:Nx-1
        xd = (x[i]+x[i+1])/2
        @inbounds for j in 1:Ny-1
            yd= (y[j]+y[j+1])/2

            B1xdy[i, j] = BxBeta(xd, y[j], beta)
            B2xdy[i, j] = ByBeta(xd, y[j], beta)
            B1xyd[i, j] = BxBeta(x[i], yd, beta)
            B2xyd[i, j] = ByBeta(x[i], yd, beta)
            B1xdyd[i, j] = BxBeta(xd, yd, beta)
            B2xdyd[i, j] = ByBeta(xd, yd, beta)
        end
    end

    for i in 1:Nx
        @inbounds for j in 1:Ny
            Numl = lexico(i, j, Nxy)
            X[i, j] = x[i]
            Y[i, j] = y[j]
            Matb1[i, j] = BxBeta(x[i], y[j], beta)
            Matb2[i, j] = ByBeta(x[i], y[j], beta)

            SolM[i, j] = SolBeta(x[i], y[j], beta, epsilon, nPeriod)
            LapParaSolM[i, j] = LapParaSolBeta(x[i], y[j], beta, epsilon, nPeriod)
            LapOrthoSolM[i,j] = LapOrthoSolBeta(x[i], y[j], beta, epsilon, nPeriod)
            SolV[Numl] = SolM[i, j]
            RhsM[i, j] = RhsBeta(x[i], y[j], beta, epsilon, nPeriod)/epsilon
            RhsV[Numl] = RhsM[i, j]
            bGradSolM[i, j] = bGradSolBeta(x[i], y[j], beta, epsilon, nPeriod)
    
        end
    end

    # Vectorisation : syntax dot
    Mat11 = Matb1.*Matb1;
    Mat12 = Matb1.*Matb2;
    Mat22 = Matb2.*Matb2;
    Mat21 = Mat12;

    Matb11xdy = B1xdy.* B1xdy
    Matb12xdy = B1xdy.* B2xdy
    Matb22xyd = B2xyd.* B2xyd
    Matb21xyd = B1xyd.* B2xyd;
    Mat21xdyd = B1xdyd .* B2xdyd;

    MatOnes = ones(Nx,Ny)
    MatZeros = zeros(Nx,Ny)

    # Flux = Sym
    LapPara = LapParaSym(Mat11,Mat12,Mat21,Mat22, Dx, Dy, Nxy)
    Lap = LapParaSym(MatOnes,MatZeros,MatZeros,MatOnes, Dx, Dy, Nxy)

    LapOrtho = Lap - LapPara
    Mat::SparseMatrixCSC{Float64} = LapOrtho + LapPara/epsilon # Typage

    @inbounds for j in 1:Ny
        Numl = lexico(1, j, Nxy)
        Mat[Numl, :] .= 0
        Mat[Numl, Numl] = 1/(Dx*Dx)
        RhsV[Numl] = SolM[1,j]*Mat[Numl, Numl]

        Numl = lexico(Nx, j, Nxy)
        Mat[Numl, :] .= 0
        Mat[Numl, Numl] = 1/(Dx*Dx)
        RhsV[Numl] = SolM[Nx, j]*Mat[Numl, Numl]
    end 

    @inbounds for i in 2:Nx-1
        xi = x[i]
        Numl = lexico(i, 1, Nxy)
        yj = (y[1]+y[2])*0.5
        RhsV[Numl] = yTraceBeta(xi, yj, beta, epsilon, nPeriod)/Dy
        Numl = lexico(i,Ny, Nxy)
        yj = (y[Ny]+y[Ny-1])*0.5
        RhsV[Numl] = -yTraceBeta(xi, yj, beta, epsilon, nPeriod)/Dy
    end

    Approx = Mat\RhsV
    ApproxV = transpose(reshape(Approx,Ny,Nx))
    Error = ApproxV - SolM

    err =  maximum(abs.(Error))
    # println("Error = ", err)
end

@benchmark main() 
# main()