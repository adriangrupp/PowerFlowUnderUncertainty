# Initialize for (sparse) non-intrusive PCE - global parameters, polynomial basis, and samples
using Random
Random.seed!(1234)
sys = setupPowerSystem()
μ, σ, w = [2.1, 3.2], [0.3, 0.4], [0.3, 0.7]
@assert sum(w) == 1 "The weights do not sum to one."
deg = 4 # tuning parameter
unc = setupUncertaintySparse(μ,σ,w,sys[:Nd],deg)
ξ = sampleFromGaussianMixture(5000,μ,σ,w)
