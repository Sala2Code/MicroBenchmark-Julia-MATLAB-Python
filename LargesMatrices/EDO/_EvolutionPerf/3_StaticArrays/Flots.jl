using DifferentialEquations

include("./Calculs.jl") 

function condition(X, t, integrator) # evenementRenverse
    value = X[2]*(1-X[2])
end

function affect!(integrator)    
    terminate!(integrator)
end

function propagation(Ydepart::SVector{2,Float64})::Tuple{Float64, Vector{Float64}}
    tmax = 15.
    intervalleT = (0.0, tmax)

    cb = ContinuousCallback(condition, affect!)
    fun = ODEFunction(champMagNormaliseRenverse)
    prob = ODEProblem(fun, Ydepart, intervalleT, callback = cb)
    sol = solve(prob, DP5(), reltol = rTol, abstol = aTol)

    return sol.t[end], sol.u[end]
end
    