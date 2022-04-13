using PowerFlowUnderUncertainty, PowerModels, Ipopt, JuMP, TimerOutputs, JLD2

"""
30 Bus net: Sparse PCE for stochastic power flow.
2 Uncertainties: Bus 8 (Load 5) and Bus 13 (Generator 6).
Take samples of power values, compute PF, perform regression for all needed variables.
"""

caseFile = "case30.m"
numSamples = 60
maxDeg = 4
nUnc = 2
postProcessing = false

println("\n\t\t===== Stochastic Power Flow: 30 Bus case, 2 Uncertainties, sparse PCE =====\n")

## Read case file, initialize network uncertainties and corresponding values
include("init_ni.jl")
network_data = readCaseFile(caseFile)
sys = parseNetworkData(network_data)
# Define uncertain buses
p = [sys[:Pd][5], sys[:Pg][6]]
q = [sys[:Qd][5], sys[:Qg][6]]
unc = initUncertainty_2u(p, q, numSamples)

## Dict of simulation results: each row of a parameter describes a bus
pfRes = Dict(:pg => Array{Float64}(undef, sys[:Ng], 0),
    :qg => Array{Float64}(undef, sys[:Ng], 0),
    :e => Array{Float64}(undef, sys[:Nbus], 0),
    :f => Array{Float64}(undef, sys[:Nbus], 0)
)

## Initialize solver
solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 2)

## Execute the model for all samples
println("Running $numSamples deterministic PF calculations (model evalutations)...")
@timeit to "Model evaluations" begin
    for x in eachrow(unc[:samples_bus]) # each row is a sample set
        network_data["load"]["5"]["pd"] = x[1]      # bus 8 / load 5
        network_data["load"]["5"]["qd"] = x[nUnc+1] # second half of matrix are q values
        network_data["gen"]["6"]["pg"] = x[2]      # bus 13 / gen 6, no q values vor generators

        res = runPfModel(network_data, solver) # pass modified network data

        pfRes[:pg] = hcat(pfRes[:pg], res[:pg])
        pfRes[:qg] = hcat(pfRes[:qg], res[:qg])
        pfRes[:e] = hcat(pfRes[:e], res[:e])
        pfRes[:f] = hcat(pfRes[:f], res[:f])
    end
end
println("Finished.")

# Perform the regression for PCE coefficients on pd, qd, e and f
println("\nCompute sparse PCE coefficients...\n")
@timeit to "PCE Regression" begin
    pce, mse = computeCoefficientsSparse(unc[:samples_unc], pfRes, unc, K=7)
end

# Get additional PCE of currents, branch flows and demands
pf_state = getGridStateNonintrusive(pce, sys, unc)
println("PCE coefficients:")
display(pf_state)

# Sample for parameter values, using their PCE representations
pf_samples = generateSamples(unc[:ξ], pf_state, sys, unc)
println("\nPCE model evaluations for $(size(unc[:ξ])[1]) samples ξ:")
display(pf_samples)
println()


### Store data for evaluation ###

# PCE coefficients
f_coeff = "coefficients/SPF_2u_sparse.jld2"
save(f_coeff, "pf_state", pf_state)
println("PCE coefficients data saved to $f_coeff.\n")

# Compute and store moments from PCE coefficients
@timeit to "Moments calculation" begin
    moments = computeMoments(pf_state, unc)
end

f_moms = "coefficients/SPF_2u_sparse_moments.jld2"
save(f_moms, "moments", moments)
println("PCE moments data saved to $f_moms.\n")

## Show timing stats
println("Timing resutls:")
show(to)


### POST PROCESSING ###
if postProcessing
    mycolor = "red"
    plotHistogram_6in9(pf_samples[:pg], "pg", "./plots/2u_sparse"; fignum=1 + 10, color=mycolor)
    plotHistogram_6in9(pf_samples[:qg], "qg", "./plots/2u_sparse"; fignum=2 + 10, color=mycolor)

    plotHistogram_9in9(pf_samples[:e][1:9, :], "e1", "./plots/2u_sparse"; fignum=3 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:e][10:18, :], "e2", "./plots/2u_sparse"; fignum=4 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:e][19:27, :], "e3", "./plots/2u_sparse"; fignum=5 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:e][28:30, :], "e4", "./plots/2u_sparse"; fignum=6 + 10, color=mycolor)

    plotHistogram_9in9(pf_samples[:f][1:9, :], "f1", "./plots/2u_sparse"; fignum=7 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:f][10:18, :], "f2", "./plots/2u_sparse"; fignum=8 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:f][19:27, :], "f3", "./plots/2u_sparse"; fignum=9 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:f][28:30, :], "f4", "./plots/2u_sparse"; fignum=10 + 10, color=mycolor)
    # Plot P & Q of uncertainties
    plotHistogram_2unc(pf_samples[:pd], pf_samples[:qd], ["pd_8", "qd_8", "pg_13", "qg_13"], "./plots/2u_sparse"; fignum=0 + 10, color=mycolor)
end
