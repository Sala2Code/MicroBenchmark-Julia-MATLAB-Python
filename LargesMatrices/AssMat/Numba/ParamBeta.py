import numpy as np
from scipy.sparse import csr_matrix
import time
from numba import jit, objmode

from Calculs import * # BxBeta & ByBeta - SolBeta - RhsBeta - bGradSolBeta - yTraceBeta
from Lap import * # LapPara & LapOrtho & LapParaSym & LapParaAsym & StackBuffer 
from lexico import * # lexico

@jit(nopython=False, forceobj=True)
def main():
    epsilon = 1e-2
    nPeriod = 1
    beta = 1
    k = 15

    Ly = 1
    Lx = 1

    Nx, Ny = k, k
    Dx = Lx/(Nx-2)
    # x=(-Dx/2:Dx:Lx+Dx/2) # Matlab : start:step:end
    x = np.arange(-Dx/2, Lx+Dx/2 + Dx, Dx) # +step pour inclure le dernier point 
    Dy = Ly/(Ny-2)
    y = np.arange(-Dy/2, Ly+Dy/2 + Dy, Dy)

    Nxy = [Nx, Ny]
    NL = Nx*Ny

    X = np.zeros((Nx, Ny))
    Y = np.zeros((Nx, Ny))
    Matb1 = np.zeros((Nx, Ny))
    Matb2 = np.zeros((Nx, Ny))
    SolM = np.zeros((Nx, Ny))

    dxSolM = np.zeros((Nx, Ny))
    dySolM = np.zeros((Nx, Ny))
    dxxSolM = np.zeros((Nx, Ny))
    dyySolM = np.zeros((Nx, Ny))
    dxySolM = np.zeros((Nx, Ny))
    RhsM = np.zeros((Nx, Ny))
    bGradSolM = np.zeros((Nx, Ny))
    SolV = np.zeros((NL, 1))
    RhsV = np.zeros((NL, 1))
    SolV1 = np.zeros((NL, 1))

    B1xdy = np.zeros((Nx, Ny))
    B2xdy = np.zeros((Nx, Ny))
    B1xyd = np.zeros((Nx, Ny))
    B2xyd = np.zeros((Nx, Ny))

    dxBx = np.zeros((Nx, Ny))
    dyBx = np.zeros((Nx, Ny))

    dxBy = np.zeros((Nx, Ny))
    dyBy = np.zeros((Nx, Ny))

    dxdSolM = np.zeros((Nx, Ny))
    dydSolM = np.zeros((Nx, Ny))
    BxdxU1d = np.zeros((Nx, Ny))
    BydyU1d = np.zeros((Nx, Ny))

    SolM1 = np.zeros((Nx, Ny))
    B1xdyd = np.zeros((Nx-1, Ny-1))
    B2xdyd = np.zeros((Nx-1, Ny-1))

    LapParaSolM = np.zeros((Nx, Ny))
    LapOrthoSolM = np.zeros((Nx, Ny))


    for i in range(Nx-1):
        xd = (x[i]+x[i+1])/2
        for j in range(Ny-1):
            yd = (y[j]+y[j+1])/2

            B1xdy[i, j] = BxBeta(xd, y[j], beta)
            B2xdy[i, j] = ByBeta(xd, y[j], beta)
            B1xyd[i, j] = BxBeta(x[i], yd, beta)
            B2xyd[i, j] = ByBeta(x[i], yd, beta)
            B1xdyd[i, j] = BxBeta(xd, yd, beta)
            B2xdyd[i, j] = ByBeta(xd, yd, beta)

    for i in range(Nx):
        for j in range(Ny):
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
    
    RhsVbis = RhsV
    Mat11 = np.array(Matb1)*np.array(Matb1)
    Mat12 = np.array(Matb1)*np.array(Matb2)
    Mat22 = np.array(Matb2)*np.array(Matb2)
    Mat21 =  np.array(Matb1)*np.array(Matb2)

    Matb11xdy = np.array(B1xdy)*np.array(B1xdy)
    Matb12xdy = np.array(B1xdy)*np.array(B2xdy)
    Matb22xyd = np.array(B2xyd)*np.array(B2xyd)
    Matb21xyd = np.array(B1xyd)*np.array(B2xyd)
    Mat21xdyd = np.array(B1xdyd)*np.array(B2xdyd)

    MatOnes = np.ones((Nx, Ny))
    MatZeros = np.zeros((Nx, Ny))

    # Flux Sym
    LapPara = LapParaSym(Mat11,Mat12,Mat21,Mat22, Dx, Dy, Nxy)
    Lap = LapParaSym(MatOnes,MatZeros,MatZeros,MatOnes, Dx, Dy, Nxy)

    LapOrtho = Lap - LapPara
    Mat = LapOrtho + LapPara/epsilon

    # * ----- SetBoundariesBeta -----
    for j in range(Ny):
        Numl = lexico(0, j, Nxy)
        Numc = lexico(1, j, Nxy)
        Mat[Numl, :] = 0
        Mat[Numl, Numl] = 1/(Dx*Dx)
        
        RhsV[Numl] = SolM[0,j]*Mat[Numl, Numl]

        # -1 sur les Nx
        Numl = lexico(Nx-1, j, Nxy)
        Numc = lexico(Nx-2, j, Nxy)

        Mat[Numl, :] = 0
        Mat[Numl, Numl] = 1/(Dx*Dx)
        RhsV[Numl] = SolM[Nx-1, j]*Mat[Numl, Numl]

    for i in range(1,Nx-1):
        xi = x[i]
        Numl = lexico(i, 0, Nxy) # -1
        yj = (y[0]+y[1])*0.5
        RhsV[Numl] = yTraceBeta(xi, yj, beta, epsilon, nPeriod)/Dy
        Numl = lexico(i,Ny-1, Nxy)
        yj = (y[Ny-1]+y[Ny-2])*0.5
        RhsV[Numl] = -yTraceBeta(xi, yj, beta, epsilon, nPeriod)/Dy
    # * ----------- end -----------
        Approx = np.linalg.lstsq(csr_matrix(Mat).todense() , RhsV)[0] # ! Non compatible !
        ApproxV = np.reshape(Approx, (Ny,Nx))
        Error = ApproxV - SolM

        err = np.max(abs(Error))
        # print('Error = {:3f}'.format( err ))

t=0
n=5
for i in range(n):
    t0 = time.time()
    main()
    t += time.time() - t0
t = t/n
print(t)