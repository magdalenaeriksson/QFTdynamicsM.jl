using QFTdynamics
using Serialization

###################################################################################################################################
# Set up
###################################################################################################################################
lines = readlines("etc/paths.txt")
datapath = chop(lines[2],head=5,tail=0)
postpath = chop(lines[3],head=5,tail=0)
plotpath = chop(lines[4],head=5,tail=0)
   

plots = Dict()
# 0)
#datatag1 = "phi4LambdaNLO_L600_ON1_RNone_Nx16_sdim2_Nsteps1000_dt100_M50_NstepsinMemory1002_testing"
#datatag0 = "phi4NinverseNLO_L100_ON8_RNone_Nx1_sdim1_Nsteps600_dt10_M100_NstepsinMemory602_Ninverse2"
#datatag0 = "phi4NinverseNLO_L100_ON8_RNone_Nx1_sdim1_Nsteps600_dt10_M100_NstepsinMemory602_Ninverse2CLEANEDREPOmultiandfast"
#datatag0 = "phi4NinverseLO_L100_ON4_RNone_n0_Nx32_sdim3_Nsteps250_dt100_M70_NstepsinMemory252_Tachyonic"
datatag0 = "phi4tachyonicNinverseNLO_L100_ON4_RMass_n0_Nx32_sdim3_Nsteps200_dt100_M70_NstepsinMemory202_Tachyonicfast"
thesolution = deserialize( datapath * "/" * datatag0 * "/SimSolution.jld"  )
plotdata(plots, thesolution, "c")  #create canvas
plotdata(plots, thesolution, "h", 1, "1/N NLO") # fill canvas
# 1)
#datatag1 = "phi4LambdaNLO_L600_ON1_RNone_Nx16_sdim2_Nsteps1000_dt100_M50_NstepsinMemory1002_testing"
#datatag1 = "phi4LambdaLO_L100_ON8_RNone_Nx1_sdim1_Nsteps600_dt10_M100_NstepsinMemory602_Ninverse1"
#datatag1= "phi4NinverseNLO_L100_ON20_RNone_Nx1_sdim1_Nsteps600_dt10_M100_NstepsinMemory602_Ninverse2again"
datatag1 = "phi4tachyonicNinverseLO_L100_ON4_RMass_n0_Nx32_sdim3_Nsteps200_dt100_M70_NstepsinMemory202_Tachyonic"
thesolution = deserialize( datapath * "/" * datatag1 * "/SimSolution.jld"  )
##plotdata(plots, thesolution, "c")  #create canvas
plotdata(plots, thesolution, "h", 2, "1/N LO") # fill canvas
## 2)
##datatag2 = "phi4LambdaNLO_L600_ON1_RNone_Nx16_sdim2_Nsteps1000_dt100_M50_NstepsinMemory400_testing"
##datatag2 = "phi4NinverseNLO_L100_ON8_RNone_Nx1_sdim1_Nsteps600_dt10_M100_NstepsinMemory602_Ninverse2"
#datatag2 = "phi4NinverseNLO_L100_ON20_RNone_Nx1_sdim1_Nsteps600_dt10_M100_NstepsinMemory602_Ninverse2againCORRECTED"
datatag2 = "phi4tachyonicLambdaLO_L100_ON4_RMass_n0_Nx32_sdim3_Nsteps200_dt100_M70_NstepsinMemory202_Tachyonic"
thesolution2 = deserialize( datapath * "/" * datatag2 * "/SimSolution.jld"  )
plotdata(plots, thesolution2, "h", 3, "Lambda LO") # fill canvas

# saveplots)
#dirname = "Ninverse_NLO_comparisiondifferentversions"
#dirname = "phi4NinverseNLO_L100_ON8_RNone_Nx1_sdim1_Nsteps600_dt10_M100_NstepsinMemory602_Ninverse2CLEANEDREPOmultiandfast"
dirname = "TachyonicPaper"
#dirname = "phi4NinverseNLO_L100_ON20_RNone_Nx1_sdim1_Nsteps600_dt10_M100_NstepsinMemory602_Ninverse2again"
#dirname = "phi4NinverseNLO_L100_ON8_RNone_Nx1_sdim1_Nsteps600_dt10_M100_NstepsinMemory602_Ninverse2replot"
mkpath(plotpath * "/" * dirname)
mkpath(plotpath * "/" * dirname * "/pdf")
saveplots(plots, plotpath * "/" * dirname)
