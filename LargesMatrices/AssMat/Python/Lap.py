import numpy as np
import scipy.sparse as sp
from lexico import *

def StBuffer(Buffer, i, j, Val): # StackBuffer
    Buffer.ij[i+1, j+1] = 1
    Buffer.Val[i+1, j+1] += Val
    return None


def MyQx(Buffer, Val, I, J):
    StBuffer(Buffer, I[0], J[0], -Val[0])
    StBuffer(Buffer, I[0], J[1], -Val[0])
    StBuffer(Buffer, I[1], J[0], Val[0])
    StBuffer(Buffer, I[1], J[1], Val[0])
    StBuffer(Buffer, I[0], J[0], -Val[1])
    StBuffer(Buffer, I[1], J[0], -Val[1])
    StBuffer(Buffer, I[0], J[1], Val[1])
    StBuffer(Buffer, I[1], J[1], Val[1])
    return None


def MyQy(Buffer,Val,I,J):
    StBuffer(Buffer,I[0],J[0],-Val[0])
    StBuffer(Buffer,I[0],J[1],-Val[0])
    StBuffer(Buffer,I[1], J[0], Val[0])
    StBuffer(Buffer,I[1], J[1], Val[0])
    StBuffer(Buffer,I[0],J[0],-Val[1])
    StBuffer(Buffer,I[1],J[0],-Val[1])
    StBuffer(Buffer,I[0],J[1], Val[1])
    StBuffer(Buffer,I[1],J[1], Val[1]);
    return None 


class ThisBuff:
    def __init__(self):
        self.ij = np.zeros((3, 3), dtype=np.int64)
        self.Val = np.zeros((3, 3), dtype=np.float64)


#Parallel_Laplacian Operator
def LapParaSym(Mat11,Mat12,Mat21,Mat22, Dx, Dy, Nxy):

    Nx=Nxy[0]
    Ny=Nxy[1]

    # 0.0625 = 1/16
    SqDx = 0.0625/(Dx*Dx)
    SqDy = 0.0625/(Dy*Dy)
    DxDy = 0.0625/(Dy*Dx)

    V=np.zeros((3,(Nx-2)*Ny*9))

    Buffer = ThisBuff()
    Val = np.zeros((2,1)) 
    iCpt = 0 

    # Interior Nodes (for I all Nodes for J)
    for i in range(1,(Nx-1)):
        for j in range(Ny):
            if j==0:
                # -(qy)_{i,j+1/2} =  -(qy)_{i+1/2,j+1/2} - (qy)_{i-1/2,j+1/2}
                for iMean in [0,1]: # 0
                    iM = i - iMean 

                    Val[0] = -( Mat21[iM+1,j+1] + Mat21[iM,j+1] + Mat21[iM+1,j]  +Mat21[iM,j] ) * DxDy
                    Val[1] = -( Mat22[iM+1,j+1] + Mat22[iM,j+1] + Mat22[iM+1,j]  +Mat22[iM,j] ) * SqDy
                    I = [1-iMean, -iMean]
                    J = [1,0]
                    MyQy(Buffer, Val, I, J)

            elif j==Ny-1:
                # (qy)_{i,j-1/2} = (qy)_{i+1/2,j-1/2} + (qy)_{i-1/2,j-1/2}
                for iMean in [0,1]: #1
                    iM = i - iMean

                    Val[0] = ( Mat21[iM+1,j-1]+Mat21[iM,j-1] + Mat21[iM+1,j] + Mat21[iM,j] ) * DxDy
                    Val[1] = ( Mat22[iM+1,j-1]+Mat22[iM,j-1] + Mat22[iM+1,j] + Mat22[iM,j] ) * SqDy
                    
                    I = [1-iMean, -iMean]
                    J = [0,-1]
                    MyQy(Buffer, Val, I, J)

            else:
                # Flux in the X Direction
                # (qx)_{i+1/2,j} =  (qx)_{i+1/2,j+1/2} + (qx)_{i+1/2,j-1/2}
                for jMean in [0,1]:
                    jM = j-jMean
                    Val[0] = ( Mat11[i+1,jM+1] + Mat11[i,jM+1] + Mat11[i+1,jM] + Mat11[i,jM] ) * SqDx
                    Val[1] = ( Mat12[i+1,jM+1] + Mat12[i,jM+1] + Mat12[i+1,jM] + Mat12[i,jM] ) * DxDy
                
                    I = [1,0]
                    J = [1-jMean, -jMean]
                    MyQx(Buffer, Val, I, J)
                    # -(qx)_{i-1/2,j} =  -(qx)_{i-1/2,j+1/2}-(qx)_{i-1/2,j-1/2}
                    Val[0] = -( Mat11[i-1,jM+1]+Mat11[i,jM+1] + Mat11[i-1,jM] + Mat11[i,jM] ) * SqDx
                    Val[1] = -( Mat12[i-1,jM+1]+Mat12[i,jM+1] + Mat12[i-1,jM] + Mat12[i,jM] ) * DxDy
                    I = [0,-1]
                    J =  [1-jMean, -jMean]
                    MyQx(Buffer, Val, I, J)


                # (qy)_{i,j+1/2} =  (qy)_{i+1/2,j+1/2} + (qy)_{i-1/2,j+1/2}
                for iMean in [0,1]:
                    iM = i - iMean 

                    Val[0] = ( Mat21[iM+1,j+1]+Mat21[iM,j+1] + Mat21[iM+1,j] + Mat21[iM,j] ) * DxDy
                    Val[1] = ( Mat22[iM+1,j+1]+Mat22[iM,j+1] + Mat22[iM+1,j] + Mat22[iM,j] ) * SqDy
                    I = [1-iMean, -iMean]
                    J = [1,0]
                    MyQy(Buffer,Val,I,J)
                   
                    # -(qy)_{i,j-1/2} = -(qy)_{i+1/2,j-1/2} - (qy)_{i-1/2,j-1/2}
                    Val[0] = -( Mat21[iM+1,j-1]+Mat21[iM,j-1] + Mat21[iM+1,j] + Mat21[iM,j] ) * DxDy
                    Val[1] = -( Mat22[iM+1,j-1]+Mat22[iM,j-1] + Mat22[iM+1,j] + Mat22[iM,j] ) * SqDy
                    I = [1-iMean, -iMean]
                    J = [0, -1]
                    MyQy(Buffer, Val, I, J)
                

            Numl = lexico(i, j, Nxy)
   
            for j0 in [-1,0,1]:
                jb = j0+1
                for i0 in [-1,0,1]:
                    ib = i0+1
                    if Buffer.ij[ib,jb] == 1:
                        # Store current line
                        Numc = lexico(i+i0,j+j0, Nxy)

                        V[0,iCpt] = Numl
                        V[1,iCpt] = Numc
                        V[2,iCpt] = Buffer.Val[ib,jb]
                        
                        iCpt += 1 
                        # Reset Buffer
                        Buffer.ij[ib, jb] = Buffer.Val[ib, jb] = 0

    return sp.csr_matrix((V[2,:], (V[0,:], V[1,:])), shape=(Nx*Ny, Nx*Ny))
    

