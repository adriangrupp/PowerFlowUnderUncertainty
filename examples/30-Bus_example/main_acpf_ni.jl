using PowerFlowUnderUncertainty, PowerModels, LinearAlgebra, Ipopt, DelimitedFiles, JLD

### 30 Bus net: Non-intrusive PCE for stochastic power flow ###
## Take samples of power values, compute PF, perform regression for all needed variables.
caseFile = "case30.m"
numSamples = 30
maxDeg = 3

# Initialize the network uncertainties and corresponding values
include("init_ni.jl")
initUncertainty_1(sys)


# Just for debugging
# set_optimizer_attribute(pf, "print_level", 2) # set verbosity of Ipopt output. Default is 5.
# pf = PowerModels.run_ac_pf(network_data, Ipopt.Optimizer)
# status = pf["termination_status"]
# status == OPTIMAL || status == LOCALLY_SOLVED ? nothing : error("Potentially no solution found: ", status)
# print_summary(pf["solution"])

## model(x). Wrapper function for NI-algo. Currently hard coded for bus 8
# Input:  x - sampled value for active power of PQ bus.
# Return: dict of PF outputs.
function model(x)
    p, q = x, x
    network_data["load"]["5"]["pd"] = p # bus 8 / load 5
    network_data["load"]["5"]["qd"] = q # bus 8 / load 5

    pf = PowerModels.run_pf(network_data, ACRPowerModel, Ipopt.Optimizer)

    # check if solution was feasible
    status = pf["termination_status"]
    status == OPTIMAL || status == LOCALLY_SOLVED ? nothing : error("Potentially no solution found: ", status)

    pg, qg = getPQResult(pf)
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
println("\nPCE model evaluations for samples ξ:")
display(pf_samples)
println()


### Store PCE coefficients ###
f_coeff = "coefficients/SPF_NI.jld"
save(f_coeff, "pf_state", pf_state)
println("PCE coefficients data saved to $f_coeff.\n")


### Plotting ###
mycolor = "red"
plotHistogram_gen30(pf_samples[:pg], "pg", "./plots/non-intrusive"; fignum = 3 + 10, color = mycolor)
plotHistogram_gen30(pf_samples[:qg], "qg", "./plots/non-intrusive"; fignum = 4 + 10, color = mycolor)
plotHistogram_bus(pf_samples[:pd], "pd", "./plots/non-intrusive"; fignum = 1 + 10, color = mycolor)
plotHistogram_bus(pf_samples[:qd], "qd", "./plots/non-intrusive"; fignum = 2 + 10, color = mycolor)
plotHistogram_nodal(pf_samples[:e], "e", "./plots/non-intrusive"; figbum = 5 + 10, color = mycolor)
plotHistogram_nodal(pf_samples[:f], "f", "./plots/non-intrusive"; figbum = 6 + 10, color = mycolor)


### POST PROCESSING ###
# maybe later