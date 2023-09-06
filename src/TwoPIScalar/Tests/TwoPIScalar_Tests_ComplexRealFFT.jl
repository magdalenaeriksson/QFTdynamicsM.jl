using QFTdynamics
using Serialization
using Parameters
using FFTW

#
# set up parameters
#
modelfile = TwoPIScalarfile()
parameters = Dict(
    # Discretisation
    "Mass"  =>  1,
    "Nx"    =>  6, 
    "sdim"  =>  2, 
    "dt"    =>  0.01, #0.05
    "Nsteps"=>  200,
    "NstepsinMemory"  =>  0, # 3 - min for LO
    "Nmeas" =>  0, # either 0 (all timesteps are measured) or >=3 (because the plot script wants to draw 3 graphs)
    #model
    "Mod"  =>  "phi4",
    "Lambda"=>  1,
    "ONgroup"=> 8,
    #initialisation
    #"init"  =>  "Pnr",
    #"init"  =>  "Thermal",
    "init"  =>  "TopHatT1",
    #"init"  =>  "Quench",
    #"init"  =>  "Gauss",
    "xi"    => 0.4,
    "eta"   => 5.0,
    "sig"   => 1.0,
    "T"     =>  5,
    "n"     =>  2,
    #Reno
    "Reno"  =>  "RNone",
    #"Reno"  =>  "RMass",
    #PertExp
    #"Pexp"  =>  "LambdaLO",
    #"Pexp"  =>  "LambdaNLO",
    #"Pexp"  =>  "NinverseLO",
    "Pexp"  =>  "NinverseNLO",
    # Numerics
    #"Num"   =>  "GPUreduced",
    "Num"   =>  "CPUfull",
    "Nchunks"   =>  1, # if 0 -> set to @Threads.nthreads
    #extra
    "tag"   =>  "testing", #dont use "_" in the tag
    "ov"    =>  true)

#
# create Problem
#
theproblem = QFTdynamicsProblem(modelfile, parameters)
thesolution, tmpdata = QFTdynamicsSolution(modelfile, theproblem)
initialize!(thesolution, tmpdata)

######################################################################################################
# 1st: make FFT of initial data -> it is real!
######################################################################################################
println("1) make FFT of initial data -> result is real because of symmetries in F!")

@unpack problem, simdata, measurearray = thesolution
@unpack model, pexp, disc, init, reno, num, simsetup = problem

println("   - F[1,1]")
@show simdata.F[1,1]; # eg. 5x5
Fx = rfft(simdata.F[1,1]) # -> 3x5 -- N/2+1, N
println("")
@show Fx;
Fx = real(Fx)
println("")
@show thesolution.problem.disc.ivol*brfft(Fx, thesolution.problem.disc.Nx)
println("   - done")

#println("   - F[2,2]")
#@show simdata.F[2,2]
#Fx = rfft(simdata.F[2,2])
#@show Fx
#Fx = real(Fx)
#@show thesolution.problem.disc.ivol*brfft(Fx, thesolution.problem.disc.Nx)


######################################################################################################
# 2st: compare real vs complex FFT
######################################################################################################
##
## evolve to some point
##
#evto = 10
#for i in 3:evto
#evolve!(thesolution, tmpdata, i)
#end
#
##
## check I functions at evto+1
##
#t=evto+1
#@unpack problem, simdata, measurearray = thesolution
#@unpack model, pexp, disc, init, reno, num, simsetup = problem
#expandSimData!(simdata) 
#
#simdata2 = deepcopy(simdata)
#
##work on this point
#tmone = 2
#tp = 1
#
## do stuff with complex FFT
#rx = bfft(simdata.r[tmone,tp] * thesign(tmone, tp, simdata.NstepsinMemory))
#Fx = bfft(simdata.F[tmone,tp])
#Sigrx = rx .* (Fx.^2 .- (1/12.) * rx.^2) 
#SigFx = Fx .* (Fx.^2 .-    0.75 * rx.^2) 
#simdata.Sigr[tp] .= disc.ivol^3 * (-model.Lambda^2*(model.ONgroup+2)/(6 *model.ONgroup^2)) * real( fft( Sigrx) )
#simdata.SigF[tp] .= disc.ivol^3 * (-model.Lambda^2*(model.ONgroup+2)/(18*model.ONgroup^2)) * real( fft( SigFx) )
#@show simdata.SigF[tp]
#
## do stuff with real FFT
#rx = rfft(simdata2.r[tmone,tp] * thesign(tmone, tp, simdata2.NstepsinMemory))
#Fx = rfft(simdata2.F[tmone,tp])
#Sigrx = rx .* (Fx.^2 .- (1/12.) * rx.^2) 
#SigFx = Fx .* (Fx.^2 .-    0.75 * rx.^2) 
#simdata2.Sigr[tp] .= disc.ivol^3 * (-model.Lambda^2*(model.ONgroup+2)/(6 *model.ONgroup^2)) * brfft( Sigrx, disc.Nx) 
#simdata2.SigF[tp] .= disc.ivol^3 * (-model.Lambda^2*(model.ONgroup+2)/(18*model.ONgroup^2)) * brfft( SigFx, disc.Nx) 
#@show simdata2.SigF[tp]
#
## compare both
#@show isapprox(simdata.SigF[tp],simdata2.SigF[tp], atol=1e-1 )
#@show isapprox(simdata.Sigr[tp],simdata2.Sigr[tp], atol=1e-1 )
#
#