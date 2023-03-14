@doc raw"""
    LinearSolution(sr, m_par, A, B; estim)

Calculate the linearized solution to the non-linear difference equations defined
by function [`Fsys()`](@ref), using Schmitt-Grohé & Uribe (JEDC 2004) style linearization
(apply the implicit function theorem to obtain linear observation and
state transition equations).

The Jacobian is calculated using the package `ForwardDiff`

# Arguments
- `sr`: steady-state structure (variable values, indexes, numerical parameters, ...)
- `A`,`B`: derivative of [`Fsys()`](@ref) with respect to arguments `X` [`B`] and
    `XPrime` [`A`]
- `m_par`: model parameters

# Returns
- `gx`,`hx`: observation equations [`gx`] and state transition equations [`hx`]
- `alarm_LinearSolution`,`nk`: `alarm_LinearSolution=true` when solving algorithm fails, `nk` number of
    predetermined variables
- `A`,`B`: first derivatives of [`Fsys()`](@ref) with respect to arguments `X` [`B`] and
    `XPrime` [`A`]
"""
function LinearSolution(sr::SteadyResults, m_par::ModelParameters, A::Array, B::Array; estim=false)
    
    ############################################################################
    # Check whether Steady state solves the difference equation
    ############################################################################
    length_X0 = sr.n_par.ntotal 
    X0 = zeros(length_X0) .+ ForwardDiff.Dual(0.0,0.0)
    F  = Fsys(X0, X0, sr.XSS, m_par, sr.n_par, sr.indexes)
   
    FR=realpart.(F)
    println(findall(abs.(FR).>0.001))
    println("Number of States and Controls")
    println(length(F))
    println("Max error on Fsys:")
    println(maximum(abs.(FR[:])))
    
    ############################################################################
    # Calculate Jacobians of the Difference equation F
    ############################################################################
    BA  = ForwardDiff.jacobian(x-> Fsys(x[1:length_X0], x[length_X0+1:end], sr.XSS, m_par, sr.n_par, sr.indexes,), zeros(2*length_X0))

    B   = BA[:,1:length_X0]
    A   = BA[:,length_X0+1:end]
    
    ######################################
    # Solve the linearized model: Policy Functions and LOMs
    ############################################################################
    gx, hx, alarm_LinearSolution, nk = SolveDiffEq(A, B, sr.n_par, estim)

    println("State Space Solution Done")

    return gx, hx, alarm_LinearSolution, nk, A, B
end

