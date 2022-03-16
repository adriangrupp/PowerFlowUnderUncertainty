using PowerFlowUnderUncertainty, PowerModels, LinearAlgebra, Ipopt, JuMP, JLD

### 30 Bus net: Non-intrusive PCE for stochastic power flow ###
## Take samples of power values, compute PF, perform regression for all needed variables.
caseFile = "case30.m"
numSamples = 30
maxDeg = 4
postProcessing = false

println("\n\t\t===== Stochastic Power Flow: 30 Bus case, 1 Uncertainty, non-intrusive PCE =====\n")

# Read case file, initialize network uncertainties and corresponding values
include("init_ni.jl")
initUncertainty_1(sys)

solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 2)

## model(x). Wrapper function for NI-algo. Currently hard coded for bus 8
# Input:  x - sampled value for active power of PQ bus.
# Return: dict of PF outputs.
function model(x)
    p, q = x, x
    network_data["load"]["5"]["pd"] = p # bus 8 / load 5
    network_data["load"]["5"]["qd"] = q # bus 8 / load 5

    pf = PowerModels.run_pf(network_data, ACRPowerModel, solver)

    # check if solution was feasible
    status = pf["termination_status"]
    status == OPTIMAL || status == LOCALLY_SOLVED ? nothing : error("Potentially no solution found: ", status)

    pg, qg = getGenPQResult(pf)
    vr, vi = getVoltageResult(pf)

    return Dict(:pg => pg,
        :qg => qg,
        :e => vr, # we use e and f as convention for real and imaginary voltage parts
        :f => vi)
end

## Simulation results: each row of a parameter describes a bus
pfRes = Dict(:pg => Array{Float64}(undef, N_gen, 0),
    :qg => Array{Float64}(undef, N_gen, 0),
    :e => Array{Float64}(undef, N_bus, 0),
    :f => Array{Float64}(undef, N_bus, 0))

# Execute the model for all samples
println("Running $numSamples deterministic PF calculations (model evalutations)...")
for x in samples
    res = model(x)
    pfRes[:pg] = hcat(pfRes[:pg], res[:pg])
    pfRes[:qg] = hcat(pfRes[:qg], res[:qg])
    pfRes[:e] = hcat(pfRes[:e], res[:e])
    pfRes[:f] = hcat(pfRes[:f], res[:f])
end

# Perform the actual regression for PCE coefficients on pd, qd, e and f
println("Compute non-intrusive PCE coefficients...\n")
pce = computeCoefficientsNI(samples, pfRes, unc)

# Get PCE of currents, branch flows and demands
pf_state = getGridStateNonintrusive(pce, sys, unc)
println("PCE coefficients:")
display(pf_state)

# Sample for parameter values, using their PCE representations
pf_samples = generateSamples(ξ, pf_state, sys, unc)
println("\nPCE model evaluations for $(length(ξ)) samples ξ:")
display(pf_samples)
println()


### Store data for evaluation ###

# PCE coefficients
f_coeff = "coefficients/SPF_NI.jld"
save(f_coeff, "pf_state", pf_state)
println("PCE coefficients data saved to $f_coeff.\n")

# Compute and store moments from PCE coefficients
moments = computeMoments(pf_state, unc)
f_moms = "coefficients/SPF_NI_moments.jld"
save(f_moms, "moments", moments)
println("PCE moments data saved to $f_moms.\n")



### POST PROCESSING ###
if postProcessing
    mycolor = "red"
    plotHistogram_6in9(pf_samples[:pg], "pg", "./plots/non-intrusive"; fignum = 1 + 10, color = mycolor)
    plotHistogram_6in9(pf_samples[:qg], "qg", "./plots/non-intrusive"; fignum = 2 + 10, color = mycolor)

    plotHistogram_9in9(pf_samples[:e][1:9, :], "e1", "./plots/non-intrusive"; fignum = 3 + 10, color = mycolor)
    plotHistogram_9in9(pf_samples[:e][10:18, :], "e2", "./plots/non-intrusive"; fignum = 4 + 10, color = mycolor)
    plotHistogram_9in9(pf_samples[:e][19:27, :], "e3", "./plots/non-intrusive"; fignum = 5 + 10, color = mycolor)
    plotHistogram_9in9(pf_samples[:e][28:30, :], "e4", "./plots/non-intrusive"; fignum = 6 + 10, color = mycolor)

    plotHistogram_9in9(pf_samples[:f][1:9, :], "f1", "./plots/non-intrusive"; fignum = 7 + 10, color = mycolor)
    plotHistogram_9in9(pf_samples[:f][10:18, :], "f2", "./plots/non-intrusive"; fignum = 8 + 10, color = mycolor)
    plotHistogram_9in9(pf_samples[:f][19:27, :], "f3", "./plots/non-intrusive"; fignum = 9 + 10, color = mycolor)
    plotHistogram_9in9(pf_samples[:f][28:30, :], "f4", "./plots/non-intrusive"; fignum = 10 + 10, color = mycolor)

    plotHistogram_unc(pf_samples[:pd][:], pf_samples[:qd][:], ["pd", "qd"], "./plots/non-intrusive"; fignum = 0 + 10, color = mycolor)  
end
