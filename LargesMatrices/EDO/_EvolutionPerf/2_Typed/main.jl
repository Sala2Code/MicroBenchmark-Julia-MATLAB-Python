using BenchmarkTools

include("./Calculs.jl") # Sol - champMagNormaliseRenverse - aBx - aBy
include("./Flots.jl") # propagation - champMagNormaliseRenverse - remontee

# * global var
# const aTol::Float64 = 1.e-10
# const rTol::Float64 = 1.e-10
# const beta::Int64 = 3

global aTol::Float64 = 1.e-10
global rTol::Float64 = 1.e-10
global beta::Int64 = 3

function main()
    NpointsX = 15
    NpointsY = 15

    # Static
    X = range(0, stop=1, length=NpointsX)
    Y = range(0, stop=1, length=NpointsY)

    # Pre allocations
    solTheorique = zeros(NpointsY, NpointsX)
    solPropagee = zeros(NpointsY, NpointsX)
    SolErrors = zeros(NpointsX, NpointsY)

    for i in NpointsX-1:-1:2
        for j in NpointsY-1:-1:2
            xIci = X[i]
            yIci = Y[j]    
            # changement i,j devient j,i
            solTheorique[j,i] = Sol(xIci, yIci);

            Point = [xIci, yIci]
            _, foot::Vector{Float64} = propagation(Point)
            solPropagee[j,i] = Sol(foot[1], foot[2])

            SolErrors[j, i] = abs(solTheorique[j,i] - solPropagee[j,i])
            # println(j,",",i," Sol Err:", SolErrors[j, i])
        end
    end 
end


# main()
@benchmark main()
# @profview main()
# @code_warntype main()
