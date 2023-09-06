export GenericModelA

# Phi4 theory
struct GenericModelA <: AbstractGenericModel
    omega::Float64
end

#export Renormalize
#function Renormalize(model::Phi4, pexp::PertExpansion, reno::Renormalization, disc::Discretization)
#    return model
#end
#
#function Renormalize(model::Phi4, pexp::PertExpansion, reno::RMass, disc::Discretization)
#    @unpack Lambda, ONgroup, Mass = model
#    model = Phi4(Lambda, ONgroup, RenoMass(model, pexp, disc))
#    return model
#end
#
#function RenoMass(model::Phi4, pexp::PertExpansion, disc::Discretization)
#    HMass2 = 0
#    for i in 1:length(disc.fftwhelper)
#        omega = sqrt(disc.fftwhelper[i].lev2 + disc.Mass^2)
#        HMass2 += (0.5/omega)*disc.fftwhelper[i].deg
#    end
#    HMass2 *= Hartreeprefactor(model, pexp) * (1/(disc.Nx^disc.sdim))
#    HMass2 += disc.Mass^2
#    return sqrt(HMass2)
#end