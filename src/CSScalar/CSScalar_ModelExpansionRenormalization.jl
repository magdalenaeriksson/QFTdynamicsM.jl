export Phi4

# Phi4 theory
abstract type CSScalarPhi4 <: CSScalarModel end
struct CSPhi4 <: CSScalarPhi4
    Lambda::Float64
    ONgroup::Int64
    Mass2::Float64 #Mass used for solving differnetial equation
end

abstract type CSPhi4Tachyonic <: CSScalarPhi4 end # initializes only unstable modes
struct CSPhi4TachyonicUNS <: CSPhi4Tachyonic 
    Lambda::Float64
    ONgroup::Int64
    Mass2::Float64 #Mass used for solving differnetial equation
end

struct CSPhi4TachyonicALL <: CSPhi4Tachyonic # initializes all modes
    Lambda::Float64
    ONgroup::Int64
    Mass2::Float64 #Mass used for solving differnetial equation
end

export Renormalize
function Renormalize(model::CSScalarPhi4, pexp::CSScalarPertExpansion, reno::CSScalarRNone, disc::CSScalarDiscretization)
    println("No renormalization!"); flush(stdout)
    return model
end

#function Renormalize(model::CSScalarPhi4, pexp::CSScalarPertExpansion, reno::CSScalarRMass, disc::CSScalarDiscretization)
#    @unpack Lambda, ONgroup, Mass2 = model
#    model = Phi4(Lambda, ONgroup, Mass2 + getct(model, pexp, disc) )
#    println("Renormalization: Mass2 for Evolution is (bare mass, smaller than reno mass, possibly negative):", model.Mass2); flush(stdout)
#    return model
#end
#
#function Renormalize(model::CSPhi4Tachyonic, pexp::CSScalarPertExpansion, reno::CSScalarRMass, disc::CSScalarDiscretization)
#    @unpack Lambda, ONgroup, Mass2 = model
#    println("The Mass2 before renormalization:", Mass2); flush(stdout)
#    model = CSPhi4Tachyonic(Lambda, ONgroup, -Mass2 + getct(model, pexp, disc) )
#    println("The Mass2 for Evolution is therefore (possibly negative):", model.Mass2); flush(stdout)
#    return model
#end

#function getct(model::CSScalarPhi4, pexp::CSScalarPertExpansion, disc::CSScalarDiscretization)
#    ct = 0
#    for i in 1:length(disc.fftwhelper)
#        omega = sqrt(disc.fftwhelper[i].lev2 + disc.Mass^2)
#        ct -= (0.5/omega)*disc.fftwhelper[i].deg
#    end
#    println("Counterterm (without lambda and N factor) " , ct*disc.ivol); flush(stdout)
#    ct *= model.Lambda * (model.ONgroup+2)/(6*model.ONgroup) * disc.ivol
#    return ct
#end