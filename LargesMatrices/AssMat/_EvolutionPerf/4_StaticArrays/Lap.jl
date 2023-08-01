using SparseArrays
using StaticArrays
include("./lexico.jl")

struct ThisBuff     
    ij :: MMatrix{3,3, Int64}
    Val :: MMatrix{3,3, Float64}
end

# Les fonctions mutables évite la dés(allocatons) des tableaux
function StBuff!(Buffer::ThisBuff,i,j,Val) # StackBuffer
    Buffer.ij[i+2,j+2] = 1
    Buffer.Val[i+2,j+2] += Val
    return nothing
end

function MyQx!(Buffer::ThisBuff,Val::MVector{2,Float64},I::SVector{2,Int64},J::SVector{2,Int64})
    # {i+1,j+1/2}
    StBuff!(Buffer,I[1],J[1],-Val[1])
    StBuff!(Buffer,I[1],J[2],-Val[1])
    # {i,j+1/2}
    StBuff!(Buffer,I[2],J[1], Val[1])
    StBuff!(Buffer,I[2],J[2], Val[1])        
    # {i+1/2,j+1}
    StBuff!(Buffer,I[1],J[1],-Val[2])
    StBuff!(Buffer,I[2],J[1],-Val[2])
    # {i+1/2,j  }
    StBuff!(Buffer,I[1],J[2], Val[2])
    StBuff!(Buffer,I[2],J[2], Val[2])
    return nothing
end
function MyQy!(Buffer::ThisBuff,Val::MVector{2,Float64},I::SVector{2,Int64},J::SVector{2,Int64})
    StBuff!(Buffer,I[1],J[1],-Val[1])
    StBuff!(Buffer,I[1],J[2],-Val[1])

    StBuff!(Buffer,I[2], J[1], Val[1])
    StBuff!(Buffer,I[2], J[2], Val[1])
            
    # {i+1/2,j+1}
    StBuff!(Buffer,I[1],J[1],-Val[2])
    StBuff!(Buffer,I[2],J[1],-Val[2])
    # {i+1/2,j  }
    StBuff!(Buffer,I[1],J[2], Val[2])
    StBuff!(Buffer,I[2],J[2], Val[2])
end


function LapParaSym(Mat11,Mat12,Mat21,Mat22, Dx, Dy, Nxy)
    Nx=Nxy[1][1]
    Ny=Nxy[2][1]


    # 1/16 = 0.0625 
    SqDx = 0.0625/(Dx*Dx)
    SqDy = 0.0625/(Dy*Dy)
    DxDy = 0.0625/(Dy*Dx)

    V=zeros(3,(Nx-2)*Ny*9)

    # Mieux qu'un dictionnaire
    Buffer = ThisBuff(zeros(Float64,3,3),zeros(Int64,3,3))

    # Static & mutable array
    Val = MVector{2, Float64}(0.,0.);

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
                    I = @SVector [1-iMean, -iMean]
                    J = @SVector [1,0]
                    MyQy!(Buffer,Val, I, J)
                end
            elseif  j==Ny
                # (qy)_{i,j-1/2} = (qy)_{i+1/2,j-1/2} + (qy)_{i-1/2,j-1/2}
                for iMean in 0:1 #1
                    iM = i - iMean
                    Val[1] = ( Mat21[iM+1,j-1]+Mat21[iM,j-1] + Mat21[iM+1,j] + Mat21[iM,j] ) * DxDy
                    Val[2] = ( Mat22[iM+1,j-1]+Mat22[iM,j-1] + Mat22[iM+1,j] + Mat22[iM,j] ) * SqDy
                    I=@SVector [1-iMean, -iMean]
                    J=@SVector [0,-1]
                    MyQy!(Buffer, Val, I, J)
                end
            else
                # Flux in the X Direction
                # (qx)_{i+1/2,j} =  (qx)_{i+1/2,j+1/2} + (qx)_{i+1/2,j-1/2}
                for jMean in 0:1
                    jM = j-jMean
                    Val[1] = ( Mat11[i+1,jM+1] + Mat11[i,jM+1] + Mat11[i+1,jM] + Mat11[i,jM] ) * SqDx
                    Val[2] = ( Mat12[i+1,jM+1] + Mat12[i,jM+1] + Mat12[i+1,jM] + Mat12[i,jM] ) * DxDy
                    I = @SVector [1,0]
                    J = @SVector [1-jMean, -jMean]
                    MyQx!(Buffer,Val, I, J)
                    # -(qx)_{i-1/2,j} =  -(qx)_{i-1/2,j+1/2}-(qx)_{i-1/2,j-1/2}
                    Val[1] = -( Mat11[i-1,jM+1]+Mat11[i,jM+1] + Mat11[i-1,jM] + Mat11[i,jM] ) * SqDx
                    Val[2] = -( Mat12[i-1,jM+1]+Mat12[i,jM+1] + Mat12[i-1,jM] + Mat12[i,jM] ) * DxDy
                    I= @SVector [0,-1]
                    J = @SVector [1-jMean, -jMean]
                    MyQx!(Buffer,Val, I, J)
                end
                # (qy)_{i,j+1/2} =  (qy)_{i+1/2,j+1/2} + (qy)_{i-1/2,j+1/2}
                for iMean in 0:1
                    iM = i - iMean 
                    Val[1] = ( Mat21[iM+1,j+1]+Mat21[iM,j+1] + Mat21[iM+1,j] + Mat21[iM,j] ) * DxDy
                    Val[2] = ( Mat22[iM+1,j+1]+Mat22[iM,j+1] + Mat22[iM+1,j] + Mat22[iM,j] ) * SqDy
                    I = @SVector [1-iMean, -iMean]
                    J = @SVector [1,0]
                    MyQy!(Buffer, Val, I, J)
                    # -(qy)_{i,j-1/2} = -(qy)_{i+1/2,j-1/2} - (qy)_{i-1/2,j-1/2}
                    Val[1] = -( Mat21[iM+1,j-1]+Mat21[iM,j-1] + Mat21[iM+1,j] + Mat21[iM,j] ) * DxDy
                    Val[2] = -( Mat22[iM+1,j-1]+Mat22[iM,j-1] + Mat22[iM+1,j] + Mat22[iM,j] ) * SqDy
                    I = @SVector [1-iMean, -iMean]
                    J = @SVector [0,-1]
                    MyQy!(Buffer, Val, I, J)
                end
            end
            Numl = lexico(i, j, Nxy)
            for j0 in -1:1
                jb = j0+2
                for i0 in -1:1
                    ib = i0+2
                    if Buffer.ij[ib,jb] == 1
                        # Store current line
                        Numc = lexico(i+i0,j+j0, Nxy)

                        V[1,iCpt] = Numl
                        V[2,iCpt] = Numc
                        V[3,iCpt] = Buffer.Val[ib,jb]
                        iCpt += 1 

                        # Reset Buffer
                        Buffer.ij[ib, jb] = 0
                        Buffer.Val[ib, jb] = 0
                    end
                end
            end
        end
    end

    return sparse(V[1,1:iCpt-1], V[2,1:iCpt-1], V[3,1:iCpt-1], Nx*Ny, Nx*Ny)
end
