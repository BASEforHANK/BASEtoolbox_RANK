@doc raw"""
    Fsys(X, XPrime, Xss, m_par, n_par, indexes, Γ, compressionIndexes, DC, IDC, DCD, IDCD)

Equilibrium error function: returns deviations from equilibrium around steady state.

Split computation into *Aggregate Part*, handled by [`Fsys_agg()`](@ref),
and *Heterogeneous Agent Part*.

# Arguments
- `X`,`XPrime`: deviations from steady state in periods t [`X`] and t+1 [`XPrime`]
- `Xss`: states and controls in steady state
- `Γ`, `DC`, `IDC`, `DCD`,`IDCD`: transformation matrices to retrieve marginal distributions [`Γ`],
    marginal value functions [`DC`,`IDC`], and the (linear) interpolant of the copula [`DCD`,`IDCD`] from deviations
- `indexes`,`compressionIndexes`: access `Xss` by variable names
    (DCT coefficients of compressed ``V_m`` and ``V_k`` in case of `compressionIndexes`)

# Example
```jldoctest
julia> # Solve for steady state, construct Γ,DC,IDC as in SGU()
julia> Fsys(zeros(ntotal),zeros(ntotal),XSS,m_par,n_par,indexes,Γ,compressionIndexes,DC,IDC)
*ntotal*-element Array{Float64,1}:
 0.0
 0.0
 ...
 0.0
```
"""
function Fsys(X::AbstractArray, XPrime::AbstractArray, XSS::Array{Float64,1}, m_par::ModelParameters,
              n_par::NumericalParameters, indexes::IndexStruct) 
              # The function call with Duals takes
              # Reserve space for error terms
    F = zeros(eltype(X),size(X))

    ############################################################################
    #            I. Read out argument values                                   #
    ############################################################################
    # rougly 10% of computing time, more if uncompress is actually called

    ############################################################################
    # I.1. Generate code that reads aggregate states/controls
    #      from steady state deviations. Equations take the form of:
    # r       = exp.(Xss[indexes.rSS] .+ X[indexes.r])
    # rPrime  = exp.(Xss[indexes.rSS] .+ XPrime[indexes.r])
    ############################################################################

    # @generate_equations(aggr_names)
    @generate_equations()

    LMULT       = exp.(XSS[indexes.LMULTSS] .+ X[indexes.LMULT])
    LMULTPrime  = exp.(XSS[indexes.LMULTSS] .+ XPrime[indexes.LMULT])
    Khh         = exp.(XSS[indexes.KhhSS] .+ X[indexes.Khh])
    KhhPrime  = exp.(XSS[indexes.KhhSS] .+ XPrime[indexes.Khh])
    Bhh       = exp.(XSS[indexes.Bhh] .+ X[indexes.Bhh])
    BhhPrime  = exp.(XSS[indexes.BhhSS] .+ XPrime[indexes.Bhh])

    ############################################################################
    #           III. Error term calculations (i.e. model starts here)          #
    ############################################################################

    ############################################################################
    #           III. 1. Aggregate Part #
    ############################################################################
    F                  = Fsys_agg(X, XPrime, XSS, m_par, n_par, indexes)

    # Replace het agent policy function by rep agent CEE
    F[indexes.LMULT]   = log.(LMULT) - log.(((C) - (N)^(1+m_par.γ) / (1+m_par.γ))^(-m_par.ξ))
    # Replace het agent capital supply curve by rep agent FOC
    F[indexes.Khh]     = log.(LMULT) - log.(m_par.β * APrime * m_par.ASHIFT * LMULTPrime * RBPrime / πPrime);
    F[indexes.K]       = log.(K) .- log.(Khh)
    # Replace het agent bond supply curve by no arbitrage condition
    F[indexes.Bhh]     = log(RBPrime / πPrime * APrime * m_par.ASHIFT) - log((qPrime + rPrime - 1.0) / q)
    F[indexes.B]       = log.(B) .- log.(Bhh)

    return F

end
