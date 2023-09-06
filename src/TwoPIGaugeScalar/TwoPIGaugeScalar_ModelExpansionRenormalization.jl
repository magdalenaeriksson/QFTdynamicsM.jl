export SUNgaugeScalar
export U1gaugeScalar

# SU(N) gauge-scalar theory
struct SUNgaugeScalar <: TwoPIGaugeScalarModel
    Lambda::Float64     # scalar self-coupling
    Mass  ::Float64     # scalar mass (used for solving differnetial equation)
    g     ::Float64     # gauge coupling
    N     ::Int64       # degree of group: SU(N)
    CF    ::Float64     # Casimir invariant of fundamental rep
end

# U(1) gauge-scalar theory
struct U1gaugeScalar <: TwoPIGaugeScalarModel
    Lambda::Float64 # scalar self-coupling
    Mass::Float64   # scalar mass (used for solving differnetial equation)
    g::Float64      # gauge coupling
end

export Renormalize
function Renormalize(model::SUNgaugeScalar, pexp::TwoPIGaugeScalarLoopEx, reno::TwoPIGaugeScalarRNone, disc::TwoPIGaugeScalarDiscretization)
    return model
end

function Renormalize(model::SUNgaugeScalar, pexp::TwoPIGaugeScalarLoopEx, reno::TwoPIGaugeScalarRMass, disc::TwoPIGaugeScalarDiscretization)
    @unpack Lambda, Mass, g, N = model
    model = SUNgaugeScalar(Lambda, RenoMass(model, pexp, disc), g, N, CF)
    return model
end

function RenoMass(model::SUNgaugeScalar, pexp::TwoPIGaugeScalarLoopEx, disc::TwoPIGaugeScalarDiscretization)
    M2 = 0
    for i in 1:length(disc.fftwhelper)
        omega = sqrt(disc.fftwhelper[i].lev2 + disc.Mass^2)
        M2 += (0.5/omega)*disc.fftwhelper[i].deg
    end
    M2 *= model.Lambda/2 * (1/(disc.Nx^disc.sdim))
    M2 += disc.Mass^2
    return sqrt(M2)
end