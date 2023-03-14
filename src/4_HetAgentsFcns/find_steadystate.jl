@doc raw"""
    find_steadystate(m_par)

Find the stationary equilibrium capital stock.

# Returns
- `KSS`: steady-state capital stock
- `n_par::NumericalParameters`,`m_par::ModelParameters`
"""
function find_steadystate(m_par)

# -------------------------------------------------------------------------------
# Find the stationary equilibrium for coarse grid
# -------------------------------------------------------------------------------

# Read out numerical parameters for starting guess solution with reduced income grid.
n_par   = NumericalParameters(m_par = m_par, ny = 4, nm = 10, nk = 10,  ϵ = 1e-6)

# Capital stock guesses
rSS   = (1.0 .- m_par.β)./m_par.β
capital_intensity(r) = ((r + m_par.δ_0) ./ m_par.α .* m_par.μ)^(1.0 ./ (m_par.α .- 1))
labor_supply(w) = ((1.0 .- m_par.τ_prog) .* m_par.τ_lev)^(1.0 ./ (m_par.γ .+ m_par.τ_prog)) .*
                    w^((1.0 .- m_par.τ_prog) ./ (m_par.γ .+ m_par.τ_prog))
KSS = capital_intensity(rSS) .* labor_supply(wage(capital_intensity(rSS), 1.0 ./ m_par.μ, 1.0, m_par) ./ m_par.μw)

# Write changed parameter values to n_par
n_par                   = NumericalParameters(m_par = m_par, naggrstates = length(state_names), naggrcontrols = length(control_names),
                                              aggr_names  = aggr_names, distr_names = distr_names)

return KSS, n_par, m_par

end

