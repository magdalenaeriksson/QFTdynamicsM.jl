################################################################################################################################################
# Figuring projections out in 3D
#
#   - Script starts with a bunch of examples 1.-6 (commented out)
#   - Then I create a matrix build from PL and PT operators + a random matrix 
#       I obtain the scalar factors by projecting with PL and PT, and then calculating the proportionality factor elementwise
#       If there where no random matrix the proportionality factor matix would have the same value everywhere.
#       The standard deviation of the proportionality matrix is a measure of the size of the random matrix constribution -> see plots  
#   - Then I do the same thing but instead of adding a real random matrix i add a random hermitian matrix
#       Thus mimiking the gauge propagator case. Otherwise everyhing as above, execpt the proportionality matrix is now complex.
#       But it seems like the imaginary contributions are very small. Therefore I just take the real part of the  proportionality matrix.
#
#   Note: 
#   1. Math isotropy means: Only one eigenvalue. Our "Isotropy": Two eigenvalues - one from PL two from PT
#   2. We could also get the scalar factors of PL and PT by calculating the eigenvalues. But then I have not a good measure of isotropy.
#
################################################################################################################################################
using FFTWhelper
using Plots
using LaTeXStrings
using LinearAlgebra
using QFTdynamics #for saveplots
using StatsBase

# define Nx
Nx=24
sdim=3

# for plotting
# read paths
lines = readlines("etc/paths.txt")
datapath = chop(lines[2],head=5,tail=0)
postpath = chop(lines[3],head=5,tail=0)
plotpath = chop(lines[4],head=5,tail=0)

# do plotting
plotsubdirectory = "3Dprojectors"
plots = Dict()

# set up fftw helper and set p2values in matrix
#fftwhelper = getfftwhelper(Nx,sdim)
#for i in 1:length(fftwhelper)
#    k2value = fftwhelper[i].lev2
#    for j in 1:fftwhelper[i].deg
#
#    end
#end 

P_L = Array{Array{ComplexF64}}(undef,Nx,Nx,Nx) #sdim
P_T = Array{Array{ComplexF64}}(undef,Nx,Nx,Nx) #sdim
for i in 1:Nx # "in k space"
    for j in 1:Nx # "in k space"
        for k in 1:Nx # "in k space"
            P_L[i,j,k] = zeros(3,3) 
            P_T[i,j,k] = zeros(3,3)
        end
    end
end

# initialize projectors
nvalues = zeros(Int64,Nx)
for i in 0:(Nx-1)
    if (i>Nx/2)
        nvalues[i+1] = -(Nx-i) # julia indexing
    else
        nvalues[i+1] = i  # julia indexing
    end
end

# momenta based on central derivative
#momentavalues = sin.(nvalues*(2*pi)/Nx) 
momentavalues = -im * (exp.(2*pi*im*nvalues/Nx) .- 1)
# mometa at index 0 and Int(Nx/2) + 1 are 0 (yes, there are 2).Set them manually (otherwise the Nx/2+1 is like 10^-16)
momentavalues[1] = 0
momentavalues[ Int(Nx/2) + 1 ] = 0

for i in 1:Nx # "in k space" # x component
    for j in 1:Nx # "in k space" # y component
        for k in 1:Nx # "in k space" # y component
            #p2 = momentavalues[i]^2 + momentavalues[j]^2 + momentavalues[k]^2
            p2 = momentavalues[i]*conj(momentavalues[i]) + momentavalues[j]*conj(momentavalues[j]) + momentavalues[k]*conj(momentavalues[k])
            if p2 == 0
                P_T[i,j,k][1,1] = 1 #xx component
                P_T[i,j,k][2,2] = 1 #yy component
                P_T[i,j,k][3,3] = 1 #zz component
                                    
                P_T[i,j,k][1,2] = 0 #xy component
                P_T[i,j,k][1,3] = 0 #xz component
                P_T[i,j,k][2,1] = 0 #yx component
                P_T[i,j,k][2,3] = 0 #yz component
                P_T[i,j,k][3,1] = 0 #zx component
                P_T[i,j,k][3,2] = 0 #zy component
        
                P_L[i,j,k][1,1] = 0 #xx component
                P_L[i,j,k][2,2] = 0 #yy component
                P_L[i,j,k][3,3] = 0 #zz component
                                    
                P_L[i,j,k][1,2] = 0 #xy component
                P_L[i,j,k][1,3] = 0 #xz component
                P_L[i,j,k][2,1] = 0 #yx component
                P_L[i,j,k][2,3] = 0 #yz component
                P_L[i,j,k][3,1] = 0 #zx component
                P_L[i,j,k][3,2] = 0 #zy component
           else
                # P_T
                P_T[i,j,k][1,1] = 1 - momentavalues[i] * conj(momentavalues[i]) / p2 #xx component
                P_T[i,j,k][2,2] = 1 - momentavalues[j] * conj(momentavalues[j]) / p2 #yy component
                P_T[i,j,k][3,3] = 1 - momentavalues[k] * conj(momentavalues[k]) / p2 #zz component

                P_T[i,j,k][1,2] =   - momentavalues[i] * conj(momentavalues[j]) / p2 #xy component
                P_T[i,j,k][1,3] =   - momentavalues[i] * conj(momentavalues[k]) / p2 #xz component
                P_T[i,j,k][2,1] =   - momentavalues[j] * conj(momentavalues[i]) / p2 #yx component
                P_T[i,j,k][2,3] =   - momentavalues[j] * conj(momentavalues[k]) / p2 #yz component
                P_T[i,j,k][3,1] =   - momentavalues[k] * conj(momentavalues[i]) / p2 #zx component
                P_T[i,j,k][3,2] =   - momentavalues[k] * conj(momentavalues[j]) / p2 #zy component

                # P_L
                P_L[i,j,k][1,1] = momentavalues[i] * conj(momentavalues[i]) / p2 #xx component
                P_L[i,j,k][2,2] = momentavalues[j] * conj(momentavalues[j]) / p2 #yy component
                P_L[i,j,k][3,3] = momentavalues[k] * conj(momentavalues[k]) / p2 #zz component

                P_L[i,j,k][1,2] = momentavalues[i] * conj(momentavalues[j]) / p2 #xy component
                P_L[i,j,k][1,3] = momentavalues[i] * conj(momentavalues[k]) / p2 #xz component
                P_L[i,j,k][2,1] = momentavalues[j] * conj(momentavalues[i]) / p2 #yx component
                P_L[i,j,k][2,3] = momentavalues[j] * conj(momentavalues[k]) / p2 #yz component
                P_L[i,j,k][3,1] = momentavalues[k] * conj(momentavalues[i]) / p2 #zx component
                P_L[i,j,k][3,2] = momentavalues[k] * conj(momentavalues[j]) / p2 #zy component
           end
        end
    end
end

# select one k point
ihat = 8
jhat = 3
khat = 6

###########################################################################################
## 1. Apply projector to momentum   
###########################################################################################
#println("Apply projector to momentum"); flush(stdout)
#amomentum = [ momentavalues[ihat], momentavalues[jhat], momentavalues[khat] ]
#@show P_T[ihat,jhat,khat] * amomentum # gives 0 - works
#@show P_L[ihat,jhat,khat] * amomentum # gives the momentum back - works
#
###########################################################################################
## 2. Check projector identities
###########################################################################################
#println("Check projector identities"); flush(stdout)
#@show P_T[ihat,jhat,khat]
#@show P_T[ihat,jhat,khat] * P_T[ihat,jhat,khat]
#
#@show P_L[ihat,jhat,khat]
#@show P_L[ihat,jhat,khat] * P_L[ihat,jhat,khat]
#
#@show P_L[ihat,jhat,khat] * P_T[ihat,jhat,khat]
#
###########################################################################################
## 3. Apply Ps to a random matrix
###########################################################################################
#A = rand(3,3)
#A_T = P_T[ihat,jhat,khat] * A
#@show P_T[ihat,jhat,khat]
#@show A_T
#@show A_T ./ P_T[ihat,jhat,khat]
#
#A_L = P_L[ihat,jhat,khat] * A
#@show P_L[ihat,jhat,khat]
#@show A_L
#@show A_L ./ P_L[ihat,jhat,khat]
#
#@show A .- A_L .- A_T
#
###########################################################################################
## 4. Apply Ps to a symmetric random matrix
###########################################################################################
#B = A + transpose(A)
#B_T = P_T[ihat,jhat,khat] * B
#@show P_T[ihat,jhat,khat]
#@show B_T
#@show B_T ./ P_T[ihat,jhat,khat]
#
#B_L = P_L[ihat,jhat,khat] * B
#@show P_L[ihat,jhat,khat]
#@show B_L
#@show B_L ./ P_L[ihat,jhat,khat]
#
#@show B .- B_L .- B_T
#
###########################################################################################
## 5. Apply Ps to a matrix build from P_L and P_T
###########################################################################################
#C = 2 .* P_L[ihat,jhat,khat] + 3 .* P_T[ihat,jhat,khat]
#C_T = P_T[ihat,jhat,khat] * C
#@show P_T[ihat,jhat,khat]
#@show C_T
#@show C_T ./ P_T[ihat,jhat,khat]
#
#C_L = P_L[ihat,jhat,khat] * C
#@show P_L[ihat,jhat,khat]
#@show C_L
#@show C_L ./ P_L[ihat,jhat,khat]
#
## ev of this
#eigvals(C) # contains 1xlongitudenal scalar value and 2xtransversal one
#
###########################################################################################
## 6. Apply Ps to a matrix build from P_L and P_T + a random matrix
###########################################################################################
#Drand =  0.5 .* rand(3,3)
#D = 2 .* P_L[ihat,jhat,khat] + 3 .* P_T[ihat,jhat,khat] + Drand
#D_T = P_T[ihat,jhat,khat] * D
#@show P_T[ihat,jhat,khat]
#@show D_T
## yes, they are almost the same 
#@show D_T ./ P_T[ihat,jhat,khat]
## quantify as a number and also the deviation in them
#isomat = D_T ./ P_T[ihat,jhat,khat] 
#@show mean(isomat)
#@show std(isomat)
#
#D_L = P_L[ihat,jhat,khat] * D
#@show P_L[ihat,jhat,khat]
#@show D_L
#@show D_L ./ P_L[ihat,jhat,khat]
## quantify as a number and also the deviation in them
#isomat = D_L ./ P_L[ihat,jhat,khat]
#@show mean(isomat)
#@show std(isomat)
#
##get Drand out
#@show Drand
#@show D - D_T - D_L # not the same! - Drand is somehow in std isomat
#
## ev of this
#eigvals(D)

##########################################################################################
# 7. Make fancy plot of isotropy: Apply Ps to a matrix build from P_L and P_T + a random matrix
##########################################################################################
nrvals = 20
rfactor = LinRange(0,2,nrvals) # factor of random matrix contribution (causes system to be not isotropic)
isoL     = zeros(nrvals)
isoL_std = zeros(nrvals)
isoT     = zeros(nrvals)
isoT_std = zeros(nrvals)

Lvalue = 2
Tvalue = 3
randmat = rand(3,3)
for (i, r) in enumerate(rfactor)
    D = Lvalue .* P_L[ihat,jhat,khat] + Tvalue .* P_T[ihat,jhat,khat] + r * randmat
    # T
    D_T = P_T[ihat,jhat,khat] * D
    # quantify as a number and also the deviation in them
    isomat = real(D_T ./ P_T[ihat,jhat,khat] )
    isoT[i] = mean(isomat)
    isoT_std[i] = std(isomat)

    # L
    D_L = P_L[ihat,jhat,khat] * D
    # quantify as a number and also the deviation in them
    isomat = real(D_L ./ P_L[ihat,jhat,khat])
    isoL[i] = mean(isomat)
    isoL_std[i] = std(isomat)
end

# make plots
titlestring = L"A = 2 P^L + 3 P^R + r R \in \mathbb{R}^{3x3}"
plots["PEigenvalues.png"] = plot(title=titlestring, xlabel="r", ylabel="EVs of A", xlim=(rfactor[1], rfactor[end]), ylim=(0,5), legend=:bottomright)
hline!(plots["PEigenvalues.png"], [Lvalue], label="L EV", color=:blue, linewidth=1, linestyle=:dash) 
hline!(plots["PEigenvalues.png"], [Tvalue], label="T EV", color=:red , linewidth=1, linestyle=:dash) 
plot!(plots["PEigenvalues.png"],  rfactor, isoL, yerr=isoL_std, label = "measured L EV", color=:blue, linewidth=2, linestyle=:solid) 
plot!(plots["PEigenvalues.png"],  rfactor, isoT, yerr=isoT_std, label = "measured T EV", color=:red , linewidth=2, linestyle=:solid) 

plots["PEigenvalues_std.png"] = plot(xlabel="r", ylabel="std of EVs of A", xlim=(rfactor[1], rfactor[end]), ylim=(0, false), legend=:bottomright)
plot!(plots["PEigenvalues_std.png"],  rfactor, isoL_std, label = "measured stddev L EV", color=:blue, linewidth=2, linestyle=:solid) 
plot!(plots["PEigenvalues_std.png"],  rfactor, isoT_std, label = "measured stddev T EV", color=:red , linewidth=2, linestyle=:solid) 

##########################################################################################
# 8. Make fancy plot of isotropy: Apply Ps to a hermitian matrix build from P_L and P_T + a random matrix
##########################################################################################
nrvals = 20
rfactor = LinRange(0,2,nrvals) # factor of random matrix contribution (causes system to be not isotropic)
isoL     = zeros(nrvals)
isoL_std = zeros(nrvals)
isoT     = zeros(nrvals)
isoT_std = zeros(nrvals)

Lvalue = 2
Tvalue = 3
randmat = rand(3,3) + im * rand(3,3)
hermrandmat = randmat + transpose(conj(randmat))
for (i, r) in enumerate(rfactor)
    D = Lvalue .* P_L[ihat,jhat,khat] + Tvalue .* P_T[ihat,jhat,khat] + r * hermrandmat # hermitian
    # T
    D_T = P_T[ihat,jhat,khat] * D # Not hermitian
    # quantify as a number and also the deviation in them
    isomat = D_T ./ P_T[ihat,jhat,khat] # that is now a complex matrix, not hermitian
    isoT[i] = mean( real.(isomat) )
    isoT_std[i] = std( real.(isomat) )

    # L
    D_L = P_L[ihat,jhat,khat] * D # Not hermitian
    # quantify as a number and also the deviation in them
    isomat = D_L ./ P_L[ihat,jhat,khat]
    isoL[i] = mean( real.(isomat) )
    isoL_std[i] = std( real.(isomat) )
end

# make plots
titlestring = L"A = 2 P^L + 3 P^R + r H \in \mathbb{C}^{3x3}, \quad A = A^{\dagger}"
plots["PEigenvalues_hermitian.png"] = plot(title=titlestring, xlabel="r", ylabel="EVs of A", xlim=(rfactor[1], rfactor[end]), ylim=(0,5), legend=:bottomright)
hline!(plots["PEigenvalues_hermitian.png"], [Lvalue], label="L EV", color=:blue, linewidth=1, linestyle=:dash) 
hline!(plots["PEigenvalues_hermitian.png"], [Tvalue], label="T EV", color=:red , linewidth=1, linestyle=:dash) 
plot!(plots["PEigenvalues_hermitian.png"],  rfactor, isoL, yerr=isoL_std, label = "measured L EV", color=:blue, linewidth=2, linestyle=:solid) 
plot!(plots["PEigenvalues_hermitian.png"],  rfactor, isoT, yerr=isoT_std, label = "measured T EV", color=:red , linewidth=2, linestyle=:solid) 

plots["PEigenvalues_hermitian_std.png"] = plot(xlabel="r", ylabel="std of EVs of A", xlim=(rfactor[1], rfactor[end]), ylim=(0, false), legend=:bottomright)
plot!(plots["PEigenvalues_hermitian_std.png"],  rfactor, isoL_std, label = "measured stddev L EV", color=:blue, linewidth=2, linestyle=:solid) 
plot!(plots["PEigenvalues_hermitian_std.png"],  rfactor, isoT_std, label = "measured stddev T EV", color=:red , linewidth=2, linestyle=:solid) 

if isinteractive() == true
    displayplots(plots)
end 
# save plots
mkpath(plotpath * "/" * plotsubdirectory)
mkpath(plotpath * "/" * plotsubdirectory * "/pdf")
saveplots(plots, plotpath * "/" * plotsubdirectory )

# add index.php and permissionssave plots
cp("etc/index.php", joinpath( plotpath * "/" * plotsubdirectory ,"index.php"), force=true )
cp("etc/index.php", joinpath( plotpath * "/" * plotsubdirectory * "/pdf" ,"index.php"), force=true )
chmod( plotpath * "/" * plotsubdirectory , 0o777, recursive=true)