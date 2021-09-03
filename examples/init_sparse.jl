using Random
Random.seed!(1234)
sys = setupPowerSystem()
μ, σ, w = [2.1, 3.2], [0.3, 0.4], [0.3, 0.7]
@assert sum(w) == 1 "The weights do not sum to one."
deg = 5 # tuning parameter for comaprison with full PCE
unc = setupUncertaintySparse(μ,σ,w,sys[:Nd],deg)
ξ = sampleFromGaussianMixture(5000,μ,σ,w)
