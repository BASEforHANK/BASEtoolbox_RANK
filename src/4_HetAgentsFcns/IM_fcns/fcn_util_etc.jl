#---------------------------------------------------------------------
# Basic Functions: Return on capital, Wages, Employment, Output
#---------------------------------------------------------------------

interest(K::Number, Z::Number,N::Number, m_par::ModelParameters) = Z .* m_par.α .* (K ./ N) .^(m_par.α - 1.0) .- m_par.δ_0
wage(K::Number, Z::Number,N::Number, m_par::ModelParameters)     = Z .* (1 - m_par.α) .* (K./N) .^m_par.α

employment(K::Number, Z::Number, m_par::ModelParameters)         = (Z .* (1.0 - m_par.α) .* (m_par.τ_lev .* (1.0 - m_par.τ_prog)).^(1.0 / (1.0 - m_par.τ_prog)) 
                                                                    .* K .^(m_par.α )).^((1.0 - m_par.τ_prog)./(m_par.γ + m_par.τ_prog + (m_par.α) .* (1 - m_par.τ_prog)))
output(K::Number, Z::Number,N::Number, m_par::ModelParameters)   = Z .* K .^(m_par.α) .* N .^(1 - m_par.α)                                                                    

# price of tradable stock in steady state
qΠSS_fnc(Y::Number,RB,m_par) = m_par.ωΠ.*(1.0 .- 1.0 ./ m_par.μ).*Y./(RB ./m_par.π .- 1 .+ m_par.ιΠ) + 1.0
# steady state payout to entrepreneurs
profitsSS_fnc(Y::Number,RB, m_par) = (1.0 - m_par.ωΠ).*(1.0 .- 1.0 ./ m_par.μ) .* Y .+ m_par.ιΠ .* (qΠSS_fnc(Y,RB,m_par) .-1.0)
# Valuation of liquid wealth (stock)
value_liquid(B,qΠ,qΠlag) = 1.0 .+ (qΠ .- qΠlag)./B