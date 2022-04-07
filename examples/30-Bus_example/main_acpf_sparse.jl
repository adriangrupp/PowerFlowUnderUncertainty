using PowerFlowUnderUncertainty, PowerModels, LinearAlgebra, Ipopt, JuMP, JLD

### 30 Bus net: Sparse PCE for stochastic power flow ###
## Take samples of power values, compute PF, perform regression for all needed variables.
caseFile = "case30.m"
numSamples = 30
maxDeg = 4
postProcessing = false

println("\n\t\t===== Stochastic Power Flow: 30 Bus case, 1 Uncertainty, sparse PCE =====\n")

# Read case file, initialize network uncertainties and corresponding values
include("init_ni.jl")
network_data = readCaseFile(caseFile)
sys = parseNetworkData(network_data)
p = sys[:Pd][5]
q = sys[:Qd][5]
unc = initUncertainty_Nu(p, q, numSamples)

## Simulation results: each row of a parameter describes a bus
pfRes = Dict(:pg => Array{Float64}(undef, sys[:Ng], 0),
    :qg => Array{Float64}(undef, sys[:Ng], 0),
    :e => Array{Float64}(undef, sys[:Nbus], 0),
    :f => Array{Float64}(undef, sys[:Nbus], 0)
)

solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 2)


# Execute the model for all samples
println("Running $numSamples deterministic PF calculations (model evalutations)...")
for x in eachrow(unc[:samples_bus]) # each row is a sample set
    network_data["load"]["5"]["pd"] = x[1] # bus 8 / load 5
    network_data["load"]["5"]["qd"] = x[2] # bus 8 / load 5

    res = runPfModel(network_data, solver) # pass modified network data

    pfRes[:pg] = hcat(pfRes[:pg], res[:pg])
    pfRes[:qg] = hcat(pfRes[:qg], res[:qg])
    pfRes[:e] = hcat(pfRes[:e], res[:e])
    pfRes[:f] = hcat(pfRes[:f], res[:f])
end

# Perform the regression for PCE coefficients on pd, qd, e and f
println("Compute sparse PCE coefficients...\n")
pce, mse = computeCoefficientsSparse(unc[:samples_unc], pfRes, unc; K=5)

# Get PCE of currents, branch flows and demands
pf_state = getGridStateNonintrusive(pce, sys, unc)
println("PCE coefficients:")
display(pf_state)

# Sample for parameter values, using their PCE representations
pf_samples = generateSamples(unc[:ξ], pf_state, sys, unc)
println("\nPCE model evaluations for $(length(unc[:ξ])) samples ξ:")
display(pf_samples)
println()


### Store data for evaluation ###

# PCE coefficients
f_coeff = "coefficients/SPF_sparse.jld"
save(f_coeff, "pf_state", pf_state)
println("PCE coefficients data saved to $f_coeff.\n")

# Compute and store moments from PCE coefficients
moments = computeMoments(pf_state, unc)
f_moms = "coefficients/SPF_sparse_moments.jld"
save(f_moms, "moments", moments)
println("PCE moments data saved to $f_moms.\n")



### POST PROCESSING ###
if postProcessing
    mycolor = "red"
    plotHistogram_6in9(pf_samples[:pg], "pg", "./plots/sparse"; fignum=1 + 10, color=mycolor)
    plotHistogram_6in9(pf_samples[:qg], "qg", "./plots/sparse"; fignum=2 + 10, color=mycolor)

    plotHistogram_9in9(pf_samples[:e][1:9, :], "e1", "./plots/sparse"; fignum=3 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:e][10:18, :], "e2", "./plots/sparse"; fignum=4 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:e][19:27, :], "e3", "./plots/sparse"; fignum=5 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:e][28:30, :], "e4", "./plots/sparse"; fignum=6 + 10, color=mycolor)

    plotHistogram_9in9(pf_samples[:f][1:9, :], "f1", "./plots/sparse"; fignum=7 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:f][10:18, :], "f2", "./plots/sparse"; fignum=8 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:f][19:27, :], "f3", "./plots/sparse"; fignum=9 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:f][28:30, :], "f4", "./plots/sparse"; fignum=10 + 10, color=mycolor)

    plotHistogram_unc(pf_samples[:pd][:], pf_samples[:qd][:], ["pd", "qd"], "./plots/sparse"; fignum=0 + 10, color=mycolor)
end
