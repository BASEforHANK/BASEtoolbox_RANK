# Set aggregate steady state variabel values
ASS       = 1.0
ZSS       = 1.0
ZISS      = 1.0
μSS       = m_par.μ
μwSS      = m_par.μw
τprogSS   = m_par.τ_prog
τlevSS    = m_par.τ_lev

σSS          = 1.0
τprog_obsSS  = 1.0
GshockSS     = 1.0
RshockSS     = 1.0
TprogshockSS = 1.0

SshockSS  = 1.0
# rSS       = 1.0 + interest(KSS,1.0 / m_par.μ, NSS, m_par)
RBSS      = m_par.RB
LPSS      = 1 + rSS - RBSS
LPXASS    = 1 + rSS - RBSS
ISS       = m_par.δ_0 * KSS

πSS       = 1.0
πwSS      = 1.0

BDSS      = eps()
BSS       = m_par.BhhBAR .* YSS

# Calculate taxes and government expenditures
TSS       = dot(1.0, taxrev) + av_tax_rateSS*((1.0 .- 1.0 ./ m_par.μw).*wSS.*NSS)

# println("BSS/YSS: ", BSS/YSS)

# println("StockshareSS: ",(qΠSS_fnc(YSS,m_par.RB,m_par) .- 1.0)/BSS)
BgovSS        = BSS .- qΠSS_fnc(YSS,m_par.RB,m_par) .+ 1.0
# println("BgovSS/YSS: ", BgovSS/YSS)
GSS           = TSS - (m_par.RB./m_par.π-1.0)*BgovSS
# println("GSS/YSS: ",GSS/YSS)

mcSS      = 1.0 ./ m_par.μ


firm_profitsSS = (1.0 - mcSS) .* YSS
# println("firm_profitsSS: ", firm_profitsSS)
qΠSS      = qΠSS_fnc(YSS,RBSS,m_par)
qΠlagSS   = qΠSS
RLSS      = m_par.RB

CSS       = (YSS - m_par.δ_0 * KSS - GSS - m_par.Rbar*BDSS)
LMULTSS = ((CSS) - (NSS)^(1+m_par.γ) / (1+m_par.γ))^(-m_par.ξ)

qSS       = 1.0
mcwSS     = 1.0 ./ m_par.μw
mcwwSS    = wSS * mcwSS
uSS       = 1.0
unionprofitsSS = (1.0 - mcwSS) .* wSS .* NSS

BYSS   = BSS / YSS
TYSS   = TSS / YSS
TlagSS = TSS

YlagSS = YSS
BgovlagSS = BgovSS
GlagSS = GSS
IlagSS = ISS
wlagSS = wSS
qlagSS = qSS
ClagSS = CSS
av_tax_ratelagSS = av_tax_rateSS
τproglagSS       = τprogSS

YgrowthSS = 1.0
BgovgrowthSS = 1.0
IgrowthSS = 1.0
wgrowthSS = 1.0
CgrowthSS = 1.0
TgrowthSS = 1.0
HtSS      = 1.0