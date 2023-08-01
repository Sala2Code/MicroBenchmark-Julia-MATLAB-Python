import numpy as np
from scipy.sparse import csr_matrix
import time
import cython

from Calculs import * # BxBeta & ByBeta - SolBeta - RhsBeta - bGradSolBeta - yTraceBeta
from Lap import * # LapPara & LapOrtho & LapParaSym & LapParaAsym & StackBuffer 
from lexico import * # lexico

@cython.boundscheck(False)
def main():
    epsilon = 1e-2
    nPeriod = 1
    beta = 1
    k = 50

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

    cdef double[:, :] X = np.zeros((Nx, Ny), dtype=np.double)
    cdef double[:, :] Y = np.zeros((Nx, Ny), dtype=np.double)
    cdef double[:, :] Matb1 = np.zeros((Nx, Ny), dtype=np.double)
    cdef double[:, :] Matb2 = np.zeros((Nx, Ny), dtype=np.double)
    cdef double[:, :] SolM = np.zeros((Nx, Ny), dtype=np.double)

    cdef double[:, :] dxSolM = np.zeros((Nx, Ny), dtype=np.double)
    cdef double[:, :] dySolM = np.zeros((Nx, Ny), dtype=np.double)
    cdef double[:, :] dxxSolM = np.zeros((Nx, Ny), dtype=np.double)
    cdef double[:, :] dyySolM = np.zeros((Nx, Ny), dtype=np.double)
    cdef double[:, :] dxySolM = np.zeros((Nx, Ny), dtype=np.double)
    cdef double[:, :] RhsM = np.zeros((Nx, Ny), dtype=np.double)
    cdef double[:, :] bGradSolM = np.zeros((Nx, Ny), dtype=np.double)
    cdef double[:, :] SolV = np.zeros((Nx * Ny, 1), dtype=np.double)
    cdef double[:, :] RhsV = np.zeros((Nx * Ny, 1), dtype=np.double)
    cdef double[:, :] SolV1 = np.zeros((Nx * Ny, 1), dtype=np.double)

    cdef double[:, :] B1xdy = np.zeros((Nx, Ny), dtype=np.double)
    cdef double[:, :] B2xdy = np.zeros((Nx, Ny), dtype=np.double)
    cdef double[:, :] B1xyd = np.zeros((Nx, Ny), dtype=np.double)
    cdef double[:, :] B2xyd = np.zeros((Nx, Ny), dtype=np.double)

    cdef double[:, :] dxBx = np.zeros((Nx, Ny), dtype=np.double)
    cdef double[:, :] dyBx = np.zeros((Nx, Ny), dtype=np.double)

    cdef double[:, :] dxBy = np.zeros((Nx, Ny), dtype=np.double)
    cdef double[:, :] dyBy = np.zeros((Nx, Ny), dtype=np.double)

    cdef double[:, :] dxdSolM = np.zeros((Nx, Ny), dtype=np.double)
    cdef double[:, :] dydSolM = np.zeros((Nx, Ny), dtype=np.double)
    cdef double[:, :] BxdxU1d = np.zeros((Nx, Ny), dtype=np.double)
    cdef double[:, :] BydyU1d = np.zeros((Nx, Ny), dtype=np.double)

    cdef double[:, :] SolM1 = np.zeros((Nx, Ny), dtype=np.double)
    cdef double[:, :] B1xdyd = np.zeros((Nx - 1, Ny - 1), dtype=np.double)
    cdef double[:, :] B2xdyd = np.zeros((Nx - 1, Ny - 1), dtype=np.double)

    cdef double[:, :] LapParaSolM = np.zeros((Nx, Ny), dtype=np.double)
    cdef double[:, :] LapOrthoSolM = np.zeros((Nx, Ny), dtype=np.double)


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

    cdef double[:, :] Mat11 = np.array(Matb1) * np.array(Matb1)
    cdef double[:, :] Mat12 = np.array(Matb1) * np.array(Matb2)
    cdef double[:, :] Mat22 = np.array(Matb2) * np.array(Matb2)
    cdef double[:, :] Mat21 = np.array(Matb1) * np.array(Matb2)

    cdef double[:, :] Matb11xdy = np.array(B1xdy) * np.array(B1xdy)
    cdef double[:, :] Matb12xdy = np.array(B1xdy) * np.array(B2xdy)
    cdef double[:, :] Matb22xyd = np.array(B2xyd) * np.array(B2xyd)
    cdef double[:, :] Matb21xyd = np.array(B1xyd) * np.array(B2xyd)
    cdef double[:, :] Mat21xdyd = np.array(B1xdyd) * np.array(B2xdyd)

    cdef double[:, :] MatOnes = np.ones((Nx, Ny), dtype=np.double)
    cdef double[:, :] MatZeros = np.zeros((Nx, Ny), dtype=np.double)

    # Flux Sym
    LapPara = LapParaSym(Mat11, Mat12, Mat21, Mat22, Dx, Dy, Nxy)
    Lap = LapParaSym(MatOnes, MatZeros, MatZeros, MatOnes, Dx, Dy, Nxy)

    LapOrtho = Lap - LapPara
    Mat = LapOrtho + LapPara / epsilon

    # * ----- SetBoundariesBeta -----
    cdef double xi, yj
    # Set boundaries for j = 0 and j = Ny-1
    for j in range(Ny):
        Numl = lexico(0, j, Nxy)
        Numc = lexico(1, j, Nxy)
        Mat[Numl, :] = 0
        Mat[Numl, Numl] = 1 / (Dx * Dx)
        RhsV[Numl] = SolM[0, j] * Mat[Numl, Numl]

        Numl = lexico(Nx - 1, j, Nxy)
        Numc = lexico(Nx - 2, j, Nxy)
        Mat[Numl, :] = 0
        Mat[Numl, Numl] = 1 / (Dx * Dx)
        RhsV[Numl] = SolM[Nx - 1, j] * Mat[Numl, Numl]

    # Set boundaries for i = 0 and i = Nx-1
    for i in range(1, Nx - 1):
        xi = x[i]
        Numl = lexico(i, 0, Nxy)
        yj = (y[0] + y[1]) * 0.5
        RhsV[Numl] = yTraceBeta(xi, yj, beta, epsilon, nPeriod) / Dy

        Numl = lexico(i, Ny - 1, Nxy)
        yj = (y[Ny - 1] + y[Ny - 2]) * 0.5
        RhsV[Numl] = -yTraceBeta(xi, yj, beta, epsilon, nPeriod) / Dy
    # * ----------- end -----------

    Approx = np.linalg.lstsq(csr_matrix(Mat).todense() , RhsV, rcond=None)[0]
    ApproxV = np.reshape(Approx, (Ny,Nx))
    Error = ApproxV - SolM

    err = np.max(abs(Error))
    # print('Error = {:3f}'.format( err ))

