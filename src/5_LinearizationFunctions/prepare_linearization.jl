
@doc raw"""
    prepare_linearization(KSS, VmSS, VkSS, distrSS, n_par, m_par)

Compute a number of equilibrium objects needed for linearization.

# Arguments
- `KSS`: steady-state capital stock
- `n_par::NumericalParameters`,`m_par::ModelParameters`

# Returns
- `XSS::Array{Float64,1}`, `XSSaggr::Array{Float64,1}`: steady state vectors produced by [`@writeXSS()`](@ref)
- `indexes`, `indexes_aggr`: `struct`s for accessing `XSS`,`XSSaggr` by variable names, produced by [`@make_fn()`](@ref),
        [`@make_fnaggr()`](@ref)
- `n_par::NumericalParameters`,`m_par::ModelParameters`
"""
function prepare_linearization(KSS, n_par, m_par)
    # if n_par.verbose
    #     println("Running reduction step on steady-state value functions to prepare linearization")
    # end
    # ------------------------------------------------------------------------------
    # STEP 1: Evaluate StE to calculate steady state variable values
    # ------------------------------------------------------------------------------
    # Calculate other equilibrium quantities
    Paux      = n_par.Π^1000                                              # Calculate ergodic ince distribution from transitions
    distr_y   = Paux[1, :] 
    NSS       = employment(KSS, 1.0./(m_par.μ*m_par.μw), m_par)
    rSS       = interest(KSS, 1.0./m_par.μ, NSS, m_par) + 1.0 
    wSS       = wage(KSS, 1.0./m_par.μ, NSS, m_par)
    YSS       = output(KSS, 1.0, NSS, m_par)                              # stationary income distribution
    
    profitsSS = profitsSS_fnc(YSS,m_par.RB,m_par)
    unionprofitsSS  = (1.0 .- 1.0/m_par.μw) .* wSS .* NSS
    LC              = 1.0./m_par.μw *wSS.*NSS  
    # think again about how to make RANK and HANK taxes equivalent (size of profits)
    pr_scale = 1.0
    taxrev          = (LC)-m_par.τ_lev.*(LC).^(1.0-m_par.τ_prog) .+ (profitsSS - m_par.τ_lev.*(profitsSS./pr_scale).^(1.0-m_par.τ_prog) .* pr_scale)
    incgrossaux     = (LC) .+ profitsSS
    av_tax_rateSS   = dot(1.0, taxrev)./(dot(1.0,incgrossaux))
    
    # ------------------------------------------------------------------------------
    # DO NOT DELETE OR EDIT NEXT LINE! This is needed for parser.
    # aggregate steady state marker
    # @include "../3_Model/input_aggregate_steady_state.jl"
    
    # write to XSS vector
    @writeXSS
    
    # produce indexes to access XSS etc.
    indexes               = produce_indexes(n_par)
    indexes_aggr          = produce_indexes_aggr(n_par)
    
    @set! n_par.ntotal    = 3 + n_par.naggr
    @set! n_par.nstates   = 2 + n_par.naggrstates 
    @set! n_par.ncontrols = 1 + n_par.naggrcontrols
    @set! n_par.LOMstate_save = zeros(n_par.nstates, n_par.nstates)
    @set! n_par.State2Control_save = zeros(n_par.ncontrols, n_par.nstates)
    @set! n_par.nstates_r   = copy(n_par.nstates)
    @set! n_par.ncontrols_r = copy(n_par.ncontrols)
    @set! n_par.ntotal_r    = copy(n_par.ntotal)
    @set! n_par.PRightStates= Diagonal(ones(n_par.nstates))
    @set! n_par.PRightAll   = Diagonal(ones(n_par.ntotal))

    if n_par.n_agg_eqn != n_par.naggr - length(n_par.distr_names)
        @warn("Inconsistency in number of aggregate variables and equations")
    end

    return XSS, XSSaggr, indexes, indexes, indexes_aggr, n_par, m_par   
end
