# Initialize for (sparse) non-intrusive PCE - global parameters, uncertainties, polynomial basis, and samples
using Random
Random.seed!(1234)
sys = setupPowerSystem()
deg = 3 # Maximum PCE degree
μ, σ, w = [2.1, 3.2], [0.3, 0.4], [0.3, 0.7]
@assert sum(w) == 1 "The weights do not sum to one."


# Initialization for single uncertainty
function initSingleUncertainty()
    global ξ = sampleFromGaussianMixture(5000,μ,σ,w)
    global unc = setupUncertaintySparse(μ,σ,w,sys[:Nd],deg)
end

# Initialization for multiple uncertainties
function initMultiUncertainty(numUnc::Int)
    samp = [sampleFromGaussianMixture(5000, μ, σ, w) for i in 1:numUnc]
    global ξ = hcat(samp...)
    global unc = setupUncertaintyMulti(μ, σ, w, sys[:Nd], deg, numUnc)
end