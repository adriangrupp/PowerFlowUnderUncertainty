using PowerFlowUnderUncertainty, PowerModels, LinearAlgebra, Ipopt, JuMP, TimerOutputs, JLD2

"""
5 Bus net: Non-intrusive PCE for stochastic optimal power flow
1 Uncertianty: Bus 2 (PQ bus)
Take samples of power values, compute OPF, perform regression for all needed variables.
"""

caseFile = "case5.m"
numSamples = 30
maxDeg = 3
numUnc = 1
postProcessing = true

println("\n\t\t===== Stochastic OPF: 5 Bus case, 1 Uncertainty, non-intrusive PCE =====\n")

# Read case file, initialize network uncertainties and corresponding values
include("init_ni.jl")
p = sys[:Pd][1]
q = sys[:Qd][1]
unc = initUncertainty_1u(p, q)

## Initialize solver
solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 5)

## Timer for profiling
to = TimerOutput()

## Execute the model for all samples
println("Running $numSamples deterministic OPF calculations (model evalutations)...")
@timeit to "Model evaluations" begin
    for x in eachrow(unc[:samples_bus]) # each row is a sample set
        network_data["load"]["1"]["pd"] = x[1] # bus 2 / load 1
        network_data["load"]["1"]["qd"] = x[2] # bus 2 / load 1

        res = runOpfModel(network_data, solver) # pass modified network data

        pfRes[:pg] = hcat(pfRes[:pg], res[:pg])
        pfRes[:qg] = hcat(pfRes[:qg], res[:qg])
        pfRes[:e] = hcat(pfRes[:e], res[:e])
        pfRes[:f] = hcat(pfRes[:f], res[:f])
    end
end
println("Finished.")

## Perform the regression for PCE coefficients of pd, qd, e and f and their mean squared error (mse)
println("Compute non-intrusive PCE coefficients...\n")
@timeit to "PCE Regression" begin
    pce, mse = computeCoefficientsNI(unc[:samples_unc], pfRes, unc)
end

## Get PCE of currents, branch flows and demands
pf_state = getGridStateNonintrusive(pce, sys, unc)
println("PCE coefficients:")
display(pf_state)

## Sample for parameter values, using their PCE representations
pf_samples = generateSamples(unc[:ξ], pf_state, sys, unc)
println("\nPCE model evaluations for $(length(unc[:ξ])) samples ξ:")
display(pf_samples)
println()


### Store data for evaluation ###

# PCE coefficients
f_coeff = "coefficients/SOPF_NI.jld2"
save(f_coeff, "pf_state", pf_state)
println("PCE coefficients data saved to $f_coeff.\n")

# Compute and store moments from PCE coefficients
@timeit to "Moments calculation" begin
    moments = computeMoments(pf_state, unc)
end

f_moms = "coefficients/SOPF_NI_moments.jld2"
save(f_moms, "moments", moments)
println("PCE moments data saved to $f_moms.\n")

## Show timing stats
println("Timing resutls:")
show(to)


### POST PROCESSING ###
if postProcessing
    mycolor = "red"
    plotHistogram_9in9(pf_samples[:pg], "pg", "./plots/SOPF-non-intrusive"; fignum=1 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:qg], "qg", "./plots/SOPF-non-intrusive"; fignum=2 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:e][1:5, :], "e", "./plots/SOPF-non-intrusive"; fignum=3 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:f][1:5, :], "f", "./plots/SOPF-non-intrusive"; fignum=7 + 10, color=mycolor)

    plotHistogram_unc(pf_samples[:pd][:], pf_samples[:qd][:], ["pd", "qd"], "./plots/SOPF-non-intrusive"; fignum=0 + 10, color=mycolor)
end
