using SparseArrays
using StaticArrays
include("./lexico.jl")


function StackBuffer(Buffer,i,j,Val)
    Buffer["ij"][i+2,j+2] = 1
    Buffer["Val"][i+2,j+2] += Val
    return Buffer
end

function Qx(Buffer,Val,I,J)
    # {i+1,j+1/2}
    Buffer = StackBuffer(Buffer,I[1],J[1],-Val[1])
    Buffer = StackBuffer(Buffer,I[1],J[2],-Val[1])
    # {i,j+1/2}
    Buffer = StackBuffer(Buffer,I[2],J[1], Val[1])
    Buffer = StackBuffer(Buffer,I[2],J[2], Val[1])        
    # {i+1/2,j+1}
    Buffer = StackBuffer(Buffer,I[1],J[1],-Val[2])
    Buffer = StackBuffer(Buffer,I[2],J[1],-Val[2])
    # {i+1/2,j  }
    Buffer = StackBuffer(Buffer,I[1],J[2], Val[2])
    return StackBuffer(Buffer,I[2],J[2], Val[2])
end

function Qy(Buffer,Val,I,J)
    # {i+1,j+1/2}
    Buffer = StackBuffer(Buffer,I[1],J[1],-Val[1])
    Buffer = StackBuffer(Buffer,I[1],J[2],-Val[1])
    # {i,j+1/2}
    Buffer = StackBuffer(Buffer,I[2], J[1], Val[1])
    Buffer = StackBuffer(Buffer,I[2], J[2], Val[1])  
    # {i+1/2,j+1}
    Buffer = StackBuffer(Buffer,I[1],J[1],-Val[2])
    Buffer = StackBuffer(Buffer,I[2],J[1],-Val[2])
    # {i+1/2,j  }
    Buffer = StackBuffer(Buffer,I[1],J[2], Val[2])
    return StackBuffer(Buffer,I[2],J[2], Val[2]);    
end


function LapParaSym(Mat11,Mat12,Mat21,Mat22, Dx, Dy, Nxy)
    Nx=Nxy[1]
    Ny=Nxy[2]


    # 1/16 = 0.0625 
    SqDx = 0.0625/(Dx*Dx)
    SqDy = 0.0625/(Dy*Dy)
    DxDy = 0.0625/(Dy*Dx)

    V=zeros(3,(Nx-2)*Ny*9)

    Buffer = Dict()
    Buffer["ij"] = zeros(3, 3)
    Buffer["Val"] = zeros(3, 3)
    
    Val = zeros(1,2)
    iCpt = 1

    # Interior Nodes (for I all Nodes for J)
    for i in 2:Nx-1
        for j in 1:Ny
            if j==1
                # -(qy)_{i,j+1/2} =  -(qy)_{i+1/2,j+1/2} - (qy)_{i-1/2,j+1/2}
                for iMean in 0:1 # 0
                    iM = i - iMean 

                    Val[1] = -( Mat21[iM+1,j+1] + Mat21[iM,j+1] + Mat21[iM+1,j]  +Mat21[iM,j] ) * DxDy
                    Val[2] = -( Mat22[iM+1,j+1] + Mat22[iM,j+1] + Mat22[iM+1,j]  +Mat22[iM,j] ) * SqDy
                    I = [1-iMean, -iMean]
                    J = [1,0]
                    Buffer = Qy(Buffer,Val, I, J)
                end
            elseif  j==Ny
                # (qy)_{i,j-1/2} = (qy)_{i+1/2,j-1/2} + (qy)_{i-1/2,j-1/2}
                for iMean in 0:1 #1
                    iM = i - iMean
                    Val[1] = ( Mat21[iM+1,j-1]+Mat21[iM,j-1] + Mat21[iM+1,j] + Mat21[iM,j] ) * DxDy
                    Val[2] = ( Mat22[iM+1,j-1]+Mat22[iM,j-1] + Mat22[iM+1,j] + Mat22[iM,j] ) * SqDy
                    I= [1-iMean, -iMean]
                    J= [0,-1]
                    Buffer = Qy(Buffer, Val, I, J)
                end
            else
                # Flux in the X Direction
                # (qx)_{i+1/2,j} =  (qx)_{i+1/2,j+1/2} + (qx)_{i+1/2,j-1/2}
                for jMean in 0:1
                    jM = j-jMean
                    Val[1] = ( Mat11[i+1,jM+1] + Mat11[i,jM+1] + Mat11[i+1,jM] + Mat11[i,jM] ) * SqDx
                    Val[2] = ( Mat12[i+1,jM+1] + Mat12[i,jM+1] + Mat12[i+1,jM] + Mat12[i,jM] ) * DxDy
                    I = [1,0]
                    J = [1-jMean, -jMean]
                    Buffer = Qx(Buffer,Val, I, J)
                    # -(qx)_{i-1/2,j} =  -(qx)_{i-1/2,j+1/2}-(qx)_{i-1/2,j-1/2}
                    Val[1] = -( Mat11[i-1,jM+1]+Mat11[i,jM+1] + Mat11[i-1,jM] + Mat11[i,jM] ) * SqDx
                    Val[2] = -( Mat12[i-1,jM+1]+Mat12[i,jM+1] + Mat12[i-1,jM] + Mat12[i,jM] ) * DxDy
                    I= [0,-1]
                    J = [1-jMean, -jMean]
                    Buffer = Qx(Buffer,Val, I, J)
                end
                # (qy)_{i,j+1/2} =  (qy)_{i+1/2,j+1/2} + (qy)_{i-1/2,j+1/2}
                for iMean in 0:1
                    iM = i - iMean 
                    Val[1] = ( Mat21[iM+1,j+1]+Mat21[iM,j+1] + Mat21[iM+1,j] + Mat21[iM,j] ) * DxDy
                    Val[2] = ( Mat22[iM+1,j+1]+Mat22[iM,j+1] + Mat22[iM+1,j] + Mat22[iM,j] ) * SqDy
                    I = [1-iMean, -iMean]
                    J = [1,0]
                    Buffer = Qy(Buffer, Val, I, J)
                    # -(qy)_{i,j-1/2} = -(qy)_{i+1/2,j-1/2} - (qy)_{i-1/2,j-1/2}
                    Val[1] = -( Mat21[iM+1,j-1]+Mat21[iM,j-1] + Mat21[iM+1,j] + Mat21[iM,j] ) * DxDy
                    Val[2] = -( Mat22[iM+1,j-1]+Mat22[iM,j-1] + Mat22[iM+1,j] + Mat22[iM,j] ) * SqDy
                    I = [1-iMean, -iMean]
                    J = [0,-1]
                    Buffer = Qy(Buffer, Val, I, J)
                end
            end
            Numl = lexico(i, j, Nxy)
            for j0 in -1:1
                jb = j0+2
                for i0 in -1:1
                    ib = i0+2
                    if Buffer["ij"][ib,jb] == 1
                        # Store current line
                        Numc = lexico(i+i0,j+j0, Nxy)

                        V[1,iCpt] = Numl
                        V[2,iCpt] = Numc
                        V[3,iCpt] = Buffer["Val"][ib,jb]
                        iCpt += 1 

                        # Reset Buffer
                        Buffer["ij"][ib, jb] = 0
                        Buffer["Val"][ib, jb] = 0
                    end
                end
            end
        end
    end

    return sparse(V[1,1:iCpt-1], V[2,1:iCpt-1], V[3,1:iCpt-1], Nx*Ny, Nx*Ny)
end
