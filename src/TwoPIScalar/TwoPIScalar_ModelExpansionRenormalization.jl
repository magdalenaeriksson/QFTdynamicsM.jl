export Phi4

# Phi4 theory
abstract type TwoPIScalarPhi4 <: TwoPIScalarModel end
struct Phi4 <: TwoPIScalarPhi4
    Lambda::Float64
    ONgroup::Int64
    Mass2::Float64 #Mass used for solving differnetial equation
end

struct Phi4Tachyonic <: TwoPIScalarPhi4
    Lambda::Float64
    ONgroup::Int64
    Mass2::Float64 #Mass used for solving differnetial equation
end

export Renormalize
function Renormalize(model::TwoPIScalarPhi4, pexp::TwoPIScalarPertExpansion, reno::TwoPIScalarRNone, disc::TwoPIScalarDiscretizationLattice)
    println("No renormalization!"); flush(stdout)
    return model
end

function Renormalize(model::Phi4, pexp::TwoPIScalarPertExpansion, reno::TwoPIScalarRMass, disc::TwoPIScalarDiscretizationLattice)
    @unpack Lambda, ONgroup, Mass2 = model
    model = Phi4(Lambda, ONgroup, Mass2 + getct(model, pexp, disc) )
    println("Renormalization: Mass2 for Evolution is (bare mass, smaller than reno mass, possibly negative):", model.Mass2); flush(stdout)
    return model
end

function Renormalize(model::Phi4Tachyonic, pexp::TwoPIScalarPertExpansion, reno::TwoPIScalarRMass, disc::TwoPIScalarDiscretizationLattice)
    @unpack Lambda, ONgroup, Mass2 = model
    println("The Mass2 before renormalization:", Mass2); flush(stdout)
    model = Phi4Tachyonic(Lambda, ONgroup, -Mass2 + getct(model, pexp, disc) )
    println("The Mass2 for Evolution is therefore (possibly negative):", model.Mass2); flush(stdout)
    return model
end

function getct(model::TwoPIScalarPhi4, pexp::TwoPIScalarPertExpansion, disc::TwoPIScalarDiscretization)
    ct = 0
    for i in 1:length(disc.fftwhelper)
        omega = sqrt(disc.fftwhelper[i].lev2 + disc.Mass^2)
        ct -= (0.5/omega)*disc.fftwhelper[i].deg
    end
    println("Counterterm (without lambda and N factor) " , ct*disc.ivol); flush(stdout)
    ct *= Hartreeprefactor(model, pexp) * disc.ivol
    return ct
end

function Renormalize(model::TwoPIScalarPhi4, pexp::TwoPIScalarPertExpansion, reno::TwoPIScalarRenormalization, disc::TwoPIScalarDiscretizationCont)
    println("No renormalization possible with Continuum version!"); flush(stdout)
    return model
end