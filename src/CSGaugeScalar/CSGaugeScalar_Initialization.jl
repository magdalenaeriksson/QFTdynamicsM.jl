using Parameters
using FFTW
using LinearAlgebra
using Random, Distributions
using DataFrames

function initPhiandPi!(model::CS_SUNgaugeScalar, simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, disc::CSGaugeScalarDiscretization, rng::Random.MersenneTwister)
    @unpack flatt, f, c, A, phik, piik, rhok, chik, phix, piix, rhox, chix, ftplan, iftplan, Nx, sdim = tmpdata
    Random.seed!(123) # Setting the seed
    d = Normal(0, 1)
    
    #"Float64" in rand() -->> dist = Normal(0, 1) 
    for a in 1:4
        phix[a] = rand( d, (Nx,Nx,Nx))#rand(rng, Float64, size(phik[a]))
        piix[a] = rand( d, (Nx,Nx,Nx))#rand(rng, Float64, size(piik[a]))
        phik[a] = fft(phix[a]) / sqrt(disc.vol)
        piik[a] = fft(piix[a]) / sqrt(disc.vol)
        ## add factors of omega(k) to phik/piik 
        for i in 1:length(disc.fftwhelper)
            omega = sqrt(disc.fftwhelper[i].lev2 + model.Mass^2)
            for j in 1:disc.fftwhelper[i].deg
                idx = disc.fftwhelper[i].ind[j]
                phik[a][idx] *= 1/sqrt(2*omega)
                piik[a][idx] *= sqrt(omega/2)
            end
        end
    end
    #
    ## construct Gauss matrix equation elements: f
    flatt[1] .= real.(piik[1]) .* real.(phik[4]) .- real.(piik[4]) .* real.(phik[1]) .+ real.(piik[3]) .* real.(phik[2]) .- real.(piik[2]) .* real.(phik[3]) .+ imag.(piik[1]) .* imag.(phik[4]) .- imag.(piik[4]) .* imag.(phik[1]) .+ imag.(piik[3]) .* imag.(phik[2]) .- imag.(piik[2]) .* imag.(phik[3])
    flatt[2] .= real.(piik[3]) .* real.(phik[1]) .- real.(piik[1]) .* real.(phik[3]) .+ real.(piik[4]) .* real.(phik[2]) .- real.(piik[2]) .* real.(phik[4]) .+ imag.(piik[3]) .* imag.(phik[1]) .- imag.(piik[1]) .* imag.(phik[3]) .+ imag.(piik[4]) .* imag.(phik[2]) .- imag.(piik[2]) .* imag.(phik[4])
    flatt[3] .= real.(piik[1]) .* real.(phik[2]) .- real.(piik[2]) .* real.(phik[1]) .+ real.(piik[4]) .* real.(phik[3]) .- real.(piik[3]) .* real.(phik[4]) .+ imag.(piik[1]) .* imag.(phik[2]) .- imag.(piik[2]) .* imag.(phik[1]) .+ imag.(piik[4]) .* imag.(phik[3]) .- imag.(piik[3]) .* imag.(phik[4])
    flatt[4] .= real.(piik[1]) .* real.(phik[1]) .+ real.(piik[2]) .* real.(phik[2]) .+ real.(piik[3]) .* real.(phik[3]) .+ real.(piik[4]) .* real.(phik[4]) .+ imag.(piik[1]) .* imag.(phik[1]) .+ imag.(piik[2]) .* imag.(phik[2]) .+ imag.(piik[3]) .* imag.(phik[3]) .+ imag.(piik[4]) .* imag.(phik[4])
    # exclude zero-modes from f lattice sum:
    flatt[1][1]=0
    flatt[2][1]=0
    flatt[3][1]=0
    flatt[4][1]=0
    #
    for a in 1:4
        f[a] = sum(flatt[a])
    end
    #
    ## construct Gauss matrix equation elements: A 
    A[1,1] =  real(phik[4][1])
    A[1,2] = -real(phik[3][1])
    A[1,3] =  real(phik[2][1])
    A[1,4] = -real(phik[1][1])
    #
    A[2,1] = -real(phik[3][1])
    A[2,2] = -real(phik[4][1])
    A[2,3] =  real(phik[1][1])
    A[2,4] =  real(phik[2][1])
    #
    A[3,1] =  real(phik[2][1])
    A[3,2] = -real(phik[1][1])
    A[3,3] = -real(phik[4][1])
    A[3,4] =  real(phik[3][1])
    #
    A[4,1] = real(phik[1][1])
    A[4,2] = real(phik[2][1])
    A[4,3] = real(phik[3][1])
    A[4,4] = real(phik[4][1])
    #
    ## construct Gauss matrix equation elements: c
    c = -A\f
    #@show A * c + f 
    for a in 1:4
        piik[a][1] = c[a] + 0*im
    end
    #
    ## check constraint in x-space
    #  
    for a in 1:4
        phix[a] = real.(ifft(phik[a])) * sqrt(disc.vol)
        piix[a] = real.(ifft(piik[a])) * sqrt(disc.vol)
    end
    #
    ## construct Phi
    Phi1x =  phix[3] - im * phix[4]
    Phi2x =  phix[1] + im * phix[2]
    Phi3x = -phix[1] + im * phix[2]
    Phi4x =  phix[3] + im * phix[4]
    #
    for idx in 1:disc.vol
        simdata.Phi[idx] = SMatrix{2,2}( Phi1x[idx], Phi3x[idx], Phi2x[idx], Phi4x[idx] ) /sqrt(2)
    end
    #@show simdata.Phi[1]
    ## Pi
    Pi1x =  piix[3] - im * piix[4]
    Pi2x =  piix[1] + im * piix[2]
    Pi3x = -piix[1] + im * piix[2]
    Pi4x =  piix[3] + im * piix[4]
    #
    for idx in 1:disc.vol
        simdata.Pi[idx] = SMatrix{2,2}( Pi1x[idx], Pi3x[idx], Pi2x[idx], Pi4x[idx] ) /sqrt(2)
    end
    #
    ## calc rho 
    pauli1 = [0.0 1.0; 1.0 0.0]
    pauli2 = [0.0 -1.0*im; 1.0*im 0.0]
    pauli3 = [1.0 0.0; 0 -1.0]
    pauli4 = [1.0*im 0.0; 0.0 1.0*im]
    for idx in 1:disc.vol
        rhox[1][idx] = -imag(tr(adjoint(simdata.Pi[idx]) * pauli1 * simdata.Phi[idx]))
        rhox[2][idx] = -imag(tr(adjoint(simdata.Pi[idx]) * pauli2 * simdata.Phi[idx]))
        rhox[3][idx] = -imag(tr(adjoint(simdata.Pi[idx]) * pauli3 * simdata.Phi[idx]))
        rhox[4][idx] = -imag(tr(adjoint(simdata.Pi[idx]) * pauli4 * simdata.Phi[idx]))
    end
    #@show (tr.( tuple(adjoint(simdata.Pi[1])) .* simdata.pauli .* tuple(im * simdata.Phi[1])))    # Paul's/Finns notation
    #@show -imag(tr.(  tuple(adjoint(simdata.Pi[1])) .* simdata.pauli .* tuple(simdata.Phi[1])))   # me == Paul's/Finns notation (OK)
    #
    rhox .*= -model.g
    #
    ## check (global) Gauss constraint
    ##### @show sum(rhox[1])
    ##### @show sum(rhox[2])
    ##### @show sum(rhox[3])
    ##### @show sum(rhox[4])
    
    return
end

function initMCPhiandPi!(model::CS_SUNgaugeScalar, simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, disc::CSGaugeScalarDiscretization)
    @unpack phix, piix, Pi, Phi, pauli = simdata
    @unpack rhox = tmpdata
    #@show phix[1][2]
    #@show piix[1][2]
    ### construct Phi
    Phi1x =  phix[3] - im * phix[4]
    Phi2x =  phix[1] + im * phix[2]
    Phi3x = -phix[1] + im * phix[2]
    Phi4x =  phix[3] + im * phix[4]
    #
    ## construct Pi
    Pi1x =  piix[3] - im * piix[4]
    Pi2x =  piix[1] + im * piix[2]
    Pi3x = -piix[1] + im * piix[2]
    Pi4x =  piix[3] + im * piix[4]
    #
    ## construct matrix Phi and Pi
    for idx in 1:disc.vol
        simdata.Phi[idx] = SMatrix{2,2}( Phi1x[idx], Phi3x[idx], Phi2x[idx], Phi4x[idx] ) /sqrt(2)^2
        simdata.Pi[idx]  = SMatrix{2,2}( Pi1x[idx], Pi3x[idx], Pi2x[idx], Pi4x[idx] ) /sqrt(2)^2
        #simdata.PhiPhi[idx] = abs(tr( adjoint(simdata.Phi[idx]) * simdata.Phi[idx] ))
    end
    #
    ## calc rho 
    #pauli4 = [1.0*im 0.0; 0.0 1.0*im]
    for idx in 1:disc.vol
        rhox[1][idx] = -imag(tr(adjoint(Pi[idx]) * pauli[1] * Phi[idx]))
        rhox[2][idx] = -imag(tr(adjoint(Pi[idx]) * pauli[2] * Phi[idx]))
        rhox[3][idx] = -imag(tr(adjoint(Pi[idx]) * pauli[3] * Phi[idx]))
        rhox[4][idx] = -imag(tr(adjoint(Pi[idx]) * im * Phi[idx]))
    end
    rhox .*= -model.g
    #
    ## check (global) Gauss constraint
    #@show sum(rhox[1])
    #@show sum(rhox[2])
    #@show sum(rhox[3])
    #@show sum(rhox[4])
end

function initE!(model::CS_SUNgaugeScalar, simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, disc::CSGaugeScalarDiscretization)
    @unpack Ga, GaEcont, GaRelNumerator, GaRel, Ga2Latt, GaRelLatt, rhok, chik, rhox, chix, k2values = tmpdata
    @unpack U, Ea, Ectd, pauli, vol, Nx = simdata
    ##
    ## compute chi
    for a in 1:3
        rhok[a] = fft(rhox[a])
        #@show rhok[a][1]               # this should be zero
        chik[a] = -rhok[a] ./ k2values  # real coefficients here --> chix = ifft(chik) completely real
        chik[a][1] = 0                  # set first entry to zero manually to avoid 0/0 = NaN element 
        chix[a] = real( ifft(chik[a]) )
    end
    ##
    ## compute Ea (Ea: lattice < directional index vector < colour index vector)
    for x in 1:Nx
        for y in 1:Nx
            for z in 1:Nx
                Ea[x,y,z] = SVector(SVector( chix[1][x,y,z]-chix[1][mod1(x+1,Nx),y,z], chix[2][x,y,z]-chix[2][mod1(x+1,Nx),y,z], chix[3][x,y,z]-chix[3][mod1(x+1,Nx),y,z] ),
                                    SVector( chix[1][x,y,z]-chix[1][x,mod1(y+1,Nx),z], chix[2][x,y,z]-chix[2][x,mod1(y+1,Nx),z], chix[3][x,y,z]-chix[3][x,mod1(y+1,Nx),z] ),
                                    SVector( chix[1][x,y,z]-chix[1][x,y,mod1(z+1,Nx)], chix[2][x,y,z]-chix[2][x,y,mod1(z+1,Nx)], chix[3][x,y,z]-chix[3][x,y,mod1(z+1,Nx)] ) )
                Ga[x,y,z] = SVector(0., 0., 0.)
                GaEcont[x,y,z] = SVector(0., 0., 0.)
                GaRelNumerator[x,y,z] = SVector(0., 0., 0.)
                GaRel[x,y,z] = SVector(0., 0., 0.)
                Ectd[x,y,z] = SVector( SMatrix{2,2}(0,0,0,0), SMatrix{2,2}(0,0,0,0), SMatrix{2,2}(0,0,0,0) )
            end
        end
    end
    ##
    ## compute Gauss constraint
    for x in 1:Nx
        for y in 1:Nx
            for z in 1:Nx
                GaEcont[x,y,z] = Ea[x,y,z][1]+Ea[x,y,z][2]+Ea[x,y,z][3] - ( Ea[mod1(x-1,Nx),y,z][1] + Ea[x,mod1(y-1,Nx),z][2] +  Ea[x,y,mod1(z-1,Nx)][3] )
                Ga[x,y,z] = [ GaEcont[x,y,z][1] + rhox[1][x,y,z], GaEcont[x,y,z][2] + rhox[2][x,y,z], GaEcont[x,y,z][3] + rhox[3][x,y,z] ]
                GaRelNumerator[x,y,z] = [ GaEcont[x,y,z][1] - rhox[1][x,y,z], GaEcont[x,y,z][2] - rhox[2][x,y,z], GaEcont[x,y,z][3] - rhox[3][x,y,z] ]
                GaRel[x,y,z] = Ga[x,y,z] ./ GaRelNumerator[x,y,z]
                Ga2Latt[x,y,z] = sum(Ga[x,y,z] .^2)
                GaRelLatt[x,y,z] = sum(GaRel[x,y,z] .^2)
            end
        end
    end
    ###### @show sum(Ga2Latt)/vol
    ###### @show sum(GaRelLatt)/vol

    #@show Ga[4]
    #@show Garel[223]
    #@show Garel[4]
    return
end

function initU!(model::CS_SUNgaugeScalar, simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, disc::CSGaugeScalarDiscretization)
    for idx in 1:disc.vol
        simdata.U[idx] = SVector( SMatrix{2,2}(1,0,0,1), SMatrix{2,2}(1,0,0,1), SMatrix{2,2}(1,0,0,1) )
    end
    return
end

function checklocalGauss(model::CS_SUNgaugeScalar, simdata::SU2HiggsSimData, tmpdata::SU2HiggsTmpData, disc::CSGaugeScalarDiscretization)
    @unpack Pi, Phi, U, Ea, Ectd, E0, Ga, pauli, Ta, Nx, vol = simdata
    @unpack rhok, chik, rhox, chix = tmpdata
    #
    ## compute rho
    for idx in 1:disc.vol
        rhox[1][idx] = -imag(tr(adjoint(Pi[idx]) * pauli[1] * Phi[idx]))
        rhox[2][idx] = -imag(tr(adjoint(Pi[idx]) * pauli[2] * Phi[idx]))
        rhox[3][idx] = -imag(tr(adjoint(Pi[idx]) * pauli[3] * Phi[idx]))
    end
    rhox .*= -model.g
    #
    ## compute E contributions
    EpartGauss1 = 0*Ga
    EpartGauss2x = 0*Ga
    EpartGauss2y = 0*Ga
    EpartGauss2z = 0*Ga
    EpartGauss = 0*Ga
    for x in 1:Nx
        for y in 1:Nx
            for z in 1:Nx
                EpartGauss1[x,y,z] = Ea[x,y,z][1]+Ea[x,y,z][2]+Ea[x,y,z][3]
                EpartGauss2x[x,y,z] = real(tr.( pauli .* tuple(adjoint(U[mod1(x-1,Nx),y,z][1]) * pauli[1] * U[mod1(x-1,Nx),y,z][1]) ) * Ea[mod1(x-1,Nx),y,z][1][1]
                                         + tr.( pauli .* tuple(adjoint(U[mod1(x-1,Nx),y,z][1]) * pauli[2] * U[mod1(x-1,Nx),y,z][1]) ) * Ea[mod1(x-1,Nx),y,z][1][2]
                                         + tr.( pauli .* tuple(adjoint(U[mod1(x-1,Nx),y,z][1]) * pauli[3] * U[mod1(x-1,Nx),y,z][1]) ) * Ea[mod1(x-1,Nx),y,z][1][3] )
                EpartGauss2y[x,y,z] = real(tr.( pauli .* tuple(adjoint(U[x,mod1(y-1,Nx),z][2]) * pauli[1] * U[x,mod1(y-1,Nx),z][2]) ) * Ea[x,mod1(y-1,Nx),z][2][1]
                                         + tr.( pauli .* tuple(adjoint(U[x,mod1(y-1,Nx),z][2]) * pauli[2] * U[x,mod1(y-1,Nx),z][2]) ) * Ea[x,mod1(y-1,Nx),z][2][2]
                                         + tr.( pauli .* tuple(adjoint(U[x,mod1(y-1,Nx),z][2]) * pauli[3] * U[x,mod1(y-1,Nx),z][2]) ) * Ea[x,mod1(y-1,Nx),z][2][3] )
                EpartGauss2z[x,y,z] = real(tr.( pauli .* tuple(adjoint(U[x,y,mod1(z-1,Nx)][3]) * pauli[1] * U[x,y,mod1(z-1,Nx)][3]) ) * Ea[x,y,mod1(z-1,Nx)][3][1]
                                         + tr.( pauli .* tuple(adjoint(U[x,y,mod1(z-1,Nx)][3]) * pauli[2] * U[x,y,mod1(z-1,Nx)][3]) ) * Ea[x,y,mod1(z-1,Nx)][3][2]
                                         + tr.( pauli .* tuple(adjoint(U[x,y,mod1(z-1,Nx)][3]) * pauli[3] * U[x,y,mod1(z-1,Nx)][3]) ) * Ea[x,y,mod1(z-1,Nx)][3][3] )
                EpartGauss[x,y,z] = EpartGauss1[x,y,z] - 0.5*(EpartGauss2x[x,y,z]+EpartGauss2y[x,y,z]+EpartGauss2z[x,y,z])
                Ga[x,y,z] = [ EpartGauss[x,y,z][1] + rhox[1][x,y,z], EpartGauss[x,y,z][2] + rhox[2][x,y,z], EpartGauss[x,y,z][3] + rhox[3][x,y,z] ]
            end
        end
    end
    return
end



#
## Particle numbers
export getparticlenr
function getparticlenr(init::CSGaugeScalarParticle, disc::CSGaugeScalarDiscretization)
    return init.n .* ones(length(disc.fftwhelper))
end

function getparticlenr(init::CSGaugeScalarThermal, disc::CSGaugeScalarDiscretization)
    return [ 1/(exp(sqrt(disc.fftwhelper[i].lev2 + disc.Mass^2)/(init.T*disc.Mass))-1) for i in 1:length(disc.fftwhelper)]
end

function getparticlenr(init::CSGaugeScalarTopHatT1, disc::CSGaugeScalarDiscretization)
    n = zeros(length(disc.fftwhelper))
    k2min = 0.68 * disc.sdim * disc.Mass^2
    idxmin = findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2min))[2]
    k2max = 2.04 * disc.sdim * disc.Mass^2
    idxmax = findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2max))[2]
    eta = 2
    for i in idxmin:idxmax
        n[i] = eta 
    end
    return n
end

function getparticlenr(init::CSGaugeScalarTopHatT2, disc::CSGaugeScalarDiscretization)
    n = zeros(length(disc.fftwhelper))
    k2min = 0 * disc.sdim * disc.Mass^2
    idxmin = findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2min))[2]
    k2max = 1.9 * disc.sdim * disc.Mass^2
    idxmax = findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2max))[2]
    eta = 1.85
    for i in idxmin:idxmax
        n[i] = eta 
    end
    return n
end

function getparticlenr(init::CSGaugeScalarTopHatT3, disc::CSGaugeScalarDiscretization)
    n = zeros(length(disc.fftwhelper))
    k2min = 2.04 * disc.sdim * disc.Mass^2
    idxmin = findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2min))[2]
    k2max = 2.72 * disc.sdim * disc.Mass^2
    idxmax = findmin(abs.([disc.fftwhelper[i].lev2 for i in 1:length(disc.fftwhelper)].-k2max))[2]
    eta = 1.6
    for i in idxmin:idxmax
        n[i] = eta 
    end
    return n
end

function initSimDataSerial!(model::CS_SUNgaugeScalar, simdata::Vector{SU2HiggsSimData}, tmpdata::Vector{SU2HiggsTmpData}, disc::CSGaugeScalarDiscretization, num::CSNumericsGaugeScalarCPU)
    println("Serial initialisation"); flush(stdout)

    # tmpdata's are all the same
    for i in 1:num.Runs
        initU!(model, simdata[i], tmpdata[1], disc)
        initMCPhiandPi!(model, simdata[i], tmpdata[1], disc)
        initE!(model, simdata[i], tmpdata[1], disc)
    end
end

function initSimDataParallel!(model::CS_SUNgaugeScalar, simdata::Vector{SU2HiggsSimData}, tmpdata::Vector{SU2HiggsTmpData}, disc::CSGaugeScalarDiscretization, num::CSNumericsGaugeScalarCPU)
    println("Parallel initialisation"); flush(stdout)

    # with parallel processes, we want induividual tmpdata's
    # each chunk is evaluated serially, tmpdata is overwritten for each i (run) within a chunk (process/thread)
    @Threads.threads for ichunk in 1:num.threads
        for i in num.threadranges[ichunk]
            initU!(model, simdata[i], tmpdata[ichunk], disc)
            initMCPhiandPi!(model, simdata[i], tmpdata[ichunk], disc)
            initE!(model, simdata[i], tmpdata[ichunk], disc)
        end
    end
end


## Initialise
export initialize!
function initialize!(thesolution::QFTdynamicsSolutionCSGaugeScalar, tmpdata::Vector{SU2HiggsTmpData})
    #@unpack problem, simdata, measurearray = thesolution
    @unpack problem, simdata, measurearray, measurearrayofruns = thesolution
    @unpack model, pexp, disc, init, reno, num, simsetup = problem

    if typeof(init) == QFTdynamics.CSGaugeScalarThermalMC
        samplefile = ""
    elseif typeof(init) == QFTdynamics.CSGaugeScalarParticleMC
        # read Scalar field pi and phi from file and set values (x space) in simadata.phix / .piix 
        samplefile = "/Users/magdalenaeriksson/code/2PIcode/data/MCSampledIC_Nx32_sdim3_Mass100_n0_Samples916_B42_ith2_test"
        for i in 1:num.Runs
            df = CSV.read(samplefile * "/Sample_" * string(i) * ".csv", DataFrame)
            for idx in 1:disc.vol
                # phi
                simdata[i].phix[1][idx] = df.phi1x[idx]
                simdata[i].phix[2][idx] = df.phi2x[idx]
                simdata[i].phix[3][idx] = df.phi3x[idx]
                simdata[i].phix[4][idx] = df.phi4x[idx]
                # pi
                simdata[i].piix[1][idx] = df.pii1x[idx]
                simdata[i].piix[2][idx] = df.pii2x[idx]
                simdata[i].piix[3][idx] = df.pii3x[idx]
                simdata[i].piix[4][idx] = df.pii4x[idx]
            end
        end
        initSimDataParallel!(model, simdata, tmpdata, disc, num) # does: initU!, initMCPhiandPi!, initE!
        # note: initMCPhiandPi! constructs simdata.Phi / Pii from the piix / phix read from the file. additionally it calcs rho  

    else # no sample file
        #
        # like above set simdata[i].phix[1-4] & .piix[1-4] and the let initSimDataParallel! do the rest
        #

        # get random number generator

        # Mersenne Twister
        ##if num.seed == 0 
        ##    rng = MersenneTwister( )
        ##else
        ##    rng = MersenneTwister( num.seed )
        ##end
        ##println("using seed: ", num.seed)
        # built in rng 
        Random.seed!(123) # Setting the seed
        d = Normal(0, 1)

        # Maybe let more threads work on it
        #@Threads.threads for ichunk in 1:num.threads
        #    for i in num.threadranges[ichunk]

        # get particle nr per mode
        n = getparticlenr(init, disc)
        # initialize 
        for i in 1:num.Runs
            # get temporary storage
            # x space -> real lattices
            phix = 0*[createlattice(disc.Nx, disc.sdim),createlattice(disc.Nx, disc.sdim),createlattice(disc.Nx, disc.sdim),createlattice(disc.Nx, disc.sdim)]
            piix = 0*[createlattice(disc.Nx, disc.sdim),createlattice(disc.Nx, disc.sdim),createlattice(disc.Nx, disc.sdim),createlattice(disc.Nx, disc.sdim)]
            # k space -> complex lattices
            phik = 0*[createclattice(disc.Nx, disc.sdim),createclattice(disc.Nx, disc.sdim),createclattice(disc.Nx, disc.sdim),createclattice(disc.Nx, disc.sdim)]
            piik = 0*[createclattice(disc.Nx, disc.sdim),createclattice(disc.Nx, disc.sdim),createclattice(disc.Nx, disc.sdim),createclattice(disc.Nx, disc.sdim)]
            flatt = 0*[createlattice(disc.Nx, disc.sdim),createlattice(disc.Nx, disc.sdim),createlattice(disc.Nx, disc.sdim),createlattice(disc.Nx, disc.sdim)]
            f = zeros(4)
            c = zeros(4)
            A = zeros(4,4)

            #"Float64" in rand() -->> dist = Normal(0, 1) 
            for a in 1:4
                phix[a] = rand( d, (disc.Nx,disc.Nx,disc.Nx))#rand(rng, Float64, size(phik[a]))
                piix[a] = rand( d, (disc.Nx,disc.Nx,disc.Nx))#rand(rng, Float64, size(piik[a]))
                phik[a] = fft(phix[a]) / sqrt(disc.vol)
                piik[a] = fft(piix[a]) / sqrt(disc.vol)
                ## add factors of omega(k) to phik/piik 
                for i in 1:length(disc.fftwhelper)
                    omega = sqrt(disc.fftwhelper[i].lev2 + model.Mass^2)
                    for j in 1:disc.fftwhelper[i].deg
                        idx = disc.fftwhelper[i].ind[j]
                        phik[a][idx] *= sqrt( (n[i] + 0.5)   /omega)
                        piik[a][idx] *= sqrt( (n[i] + 0.5) *omega )
                    end
                end
            end

            # construct Gauss matrix equation elements: f
            flatt[1] .= real.(piik[1]) .* real.(phik[4]) .- real.(piik[4]) .* real.(phik[1]) .+ real.(piik[3]) .* real.(phik[2]) .- real.(piik[2]) .* real.(phik[3]) .+ imag.(piik[1]) .* imag.(phik[4]) .- imag.(piik[4]) .* imag.(phik[1]) .+ imag.(piik[3]) .* imag.(phik[2]) .- imag.(piik[2]) .* imag.(phik[3])
            flatt[2] .= real.(piik[3]) .* real.(phik[1]) .- real.(piik[1]) .* real.(phik[3]) .+ real.(piik[4]) .* real.(phik[2]) .- real.(piik[2]) .* real.(phik[4]) .+ imag.(piik[3]) .* imag.(phik[1]) .- imag.(piik[1]) .* imag.(phik[3]) .+ imag.(piik[4]) .* imag.(phik[2]) .- imag.(piik[2]) .* imag.(phik[4])
            flatt[3] .= real.(piik[1]) .* real.(phik[2]) .- real.(piik[2]) .* real.(phik[1]) .+ real.(piik[4]) .* real.(phik[3]) .- real.(piik[3]) .* real.(phik[4]) .+ imag.(piik[1]) .* imag.(phik[2]) .- imag.(piik[2]) .* imag.(phik[1]) .+ imag.(piik[4]) .* imag.(phik[3]) .- imag.(piik[3]) .* imag.(phik[4])
            flatt[4] .= real.(piik[1]) .* real.(phik[1]) .+ real.(piik[2]) .* real.(phik[2]) .+ real.(piik[3]) .* real.(phik[3]) .+ real.(piik[4]) .* real.(phik[4]) .+ imag.(piik[1]) .* imag.(phik[1]) .+ imag.(piik[2]) .* imag.(phik[2]) .+ imag.(piik[3]) .* imag.(phik[3]) .+ imag.(piik[4]) .* imag.(phik[4])
            # exclude zero-modes from f lattice sum:
            flatt[1][1]=0
            flatt[2][1]=0
            flatt[3][1]=0
            flatt[4][1]=0
            # set f
            for a in 1:4
                f[a] = sum(flatt[a])
            end
            # set A
            A[1,1] =  real(phik[4][1])
            A[1,2] = -real(phik[3][1])
            A[1,3] =  real(phik[2][1])
            A[1,4] = -real(phik[1][1])

            A[2,1] = -real(phik[3][1])
            A[2,2] = -real(phik[4][1])
            A[2,3] =  real(phik[1][1])
            A[2,4] =  real(phik[2][1])

            A[3,1] =  real(phik[2][1])
            A[3,2] = -real(phik[1][1])
            A[3,3] = -real(phik[4][1])
            A[3,4] =  real(phik[3][1])

            A[4,1] = real(phik[1][1])
            A[4,2] = real(phik[2][1])
            A[4,3] = real(phik[3][1])
            A[4,4] = real(phik[4][1])

            # construct Gauss matrix equation elements: c
            c = -A\f
            #@show A * c + f 
            for a in 1:4
                piik[a][1] = c[a] + 0*im
            end
            #
            # set values in simdata
            #  
            # phi
            simdata[i].phix[1] = real.(ifft(phik[1])) * sqrt(disc.vol) #df.phi1x[idx]
            simdata[i].phix[2] = real.(ifft(phik[2])) * sqrt(disc.vol) #df.phi2x[idx]
            simdata[i].phix[3] = real.(ifft(phik[3])) * sqrt(disc.vol) #df.phi3x[idx]
            simdata[i].phix[4] = real.(ifft(phik[4])) * sqrt(disc.vol) #df.phi4x[idx]
            # pi
            simdata[i].piix[1] = real.(ifft(piik[1])) * sqrt(disc.vol) #df.pii1x[idx]
            simdata[i].piix[2] = real.(ifft(piik[2])) * sqrt(disc.vol) #df.pii2x[idx]
            simdata[i].piix[3] = real.(ifft(piik[3])) * sqrt(disc.vol) #df.pii3x[idx]
            simdata[i].piix[4] = real.(ifft(piik[4])) * sqrt(disc.vol) #df.pii4x[idx]
        end
    
        initSimDataParallel!(model, simdata, tmpdata, disc, num) # does: initU!, initMCPhiandPi!, initE!
        # note: initMCPhiandPi! constructs simdata.Phi / Pii from the piix / phix read from the file. additionally it calcs rho  
    end

    simsetup.lastEvolStep = 2
    measure!(thesolution,tmpdata,2)
end

export initializeSingle!
function initializeSingle!(thesolution::QFTdynamicsSolutionCSGaugeScalar, tmpdata::CSGaugeScalarTmpData)
    #@unpack problem, simdata, measurearray = thesolution
    @unpack problem, simdata, measurearray, measurearrayofruns = thesolution
    @unpack model, pexp, disc, init, reno, num, simsetup = problem


    if num.seed == 0 
        rng = MersenneTwister( )
    else
        rng = MersenneTwister( num.seed )
    end
    println("using seed: ", num.seed)

    initU!(model, simdata, tmpdata, disc)
    #initPhiandPi!(model, simdata, tmpdata, disc, rng)
    initMCPhiandPi!(model, simdata, tmpdata, disc)
    
    initE!(model, simdata, tmpdata, disc)
    #checklocalGauss(model, simdata, tmpdata, disc)

    simsetup.lastEvolStep = 2
    measure!(thesolution,tmpdata,2)
end