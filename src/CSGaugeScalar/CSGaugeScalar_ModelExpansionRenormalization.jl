export CS_SUNgaugeScalar
export CS_U1gaugeScalar

# SU(N) gauge-scalar theory
struct CS_SUNgaugeScalar <: CSGaugeScalarModel
    Lambda::Float64     # scalar self-coupling
    Mass  ::Float64     # scalar mass (used for solving differnetial equation)
    g     ::Float64     # gauge coupling
    N     ::Int64       # degree of group: SU(N)
    #betaGt::Float64
    #betaHt::Float64
    #CF    ::Float64     # Casimir invariant of fundamental rep
end

# U(1) gauge-scalar theory
struct CS_U1gaugeScalar <: CSGaugeScalarModel
    Lambda::Float64 # scalar self-coupling
    Mass::Float64   # scalar mass (used for solving differnetial equation)
    g::Float64      # gauge coupling
end

export Renormalize
function Renormalize(model::CS_SUNgaugeScalar, pexp::CSGaugeScalarPertExpansion, reno::CSGaugeScalarRNone, disc::CSGaugeScalarDiscretization)
    println("CS sim - no renormalisation!"); flush(stdout)
    return model
end
