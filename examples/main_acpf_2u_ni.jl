using PowerFlowUnderUncertainty, LinearAlgebra, JuMP, Ipopt, DelimitedFiles, JLD
include("powersystem.jl")
include("init_ni.jl")

### Non-intrusive PCE for stochastic power flow with 2 uncertainties ###
## Take samples of power values, compute PF, perform regression for all needed variables.

numSamples = 20
maxDegree = deg
numUnc = 2
initMultiUncertainty(numUnc)

println("Setting up PF model.")
pf = Model(with_optimizer(Ipopt.Optimizer))
set_optimizer_attribute(pf, "print_level", 2) # set verbosity of Ipopt output. Default is 5.
addCoreDeterministic!(pf, sys)
addPVBusDeterministic!(pf)

## Model: Wrapper function for NI-algo.
# Input:  x - vector of sampled values for active power of PQ buses.
# Return: dict of PF outputs.
function model(x::AbstractVector)
    p1, q1 = x[1], x[1] * 0.85
    p2, q2 = x[2], x[2] * 0.85
    sys[:P][1] = p1
    sys[:Q][1] = q1
    sys[:P][2] = p2
    sys[:Q][2] = q2
    resetPowerFlowConstraint!(pf, sys) # initialize pf model with new sample
    optimize!(pf) # actual model function

    return Dict(:pg => value.(pf[:pg]),
        :qg => value.(pf[:qg]),
        :e => value.(pf[:e]),
        :f => value.(pf[:f]))
end

## Perform full non-intrusive computations
busRes = Dict(:pg => Array{Float64}(undef, 2, 0),
    :qg => Array{Float64}(undef, 2, 0),
    :e => Array{Float64}(undef, 4, 0),
    :f => Array{Float64}(undef, 4, 0))
samples = [sampleFromGaussianMixture(numSamples, μ, σ, w) for i in 1:numUnc]
X = hcat(samples...)

# Y = model.(X)
println("Running $numSamples deterministic PF calculations (model evalutations)...")
for x in eachrow(X)
    res = model(x)
    busRes[:pg] = hcat(busRes[:pg], res[:pg])
    busRes[:qg] = hcat(busRes[:qg], res[:qg])
    busRes[:e] = hcat(busRes[:e], res[:e])
    busRes[:f] = hcat(busRes[:f], res[:f])
end

# Perform the actual regression for PCE coefficients on pd, qd, e and f
println("Compute non-intrusive PCE coefficients...\n")
pce = computeCoefficientsNI(X, busRes, unc)

# Get PCE of currents, branch flows and demands
pf_state = getGridStateNonintrusive(pce, pf, sys, unc)
println("PCE coefficients:")
display(pf_state)

# Sample for parameter values, using their PCE representations
pf_samples = generateSamples(ξ, pf_state, sys, unc)
println("\nPCE model evaluations for samples ξ:")
display(pf_samples)
println()


### Store PCE coefficients ###
f_coeff = "coefficients/SPF_NI_2unc.jld"
save(f_coeff, "pf_state", pf_state)
println("PCE coefficients data saved to $f_coeff.\n")


### Plotting ###
buscolor = "red"
nodecolor = "blue"
plotHistogram_bus(pf_samples[:pd], "pd", "./plots/non-intrusive_2unc"; fignum = 1 + 10, color = buscolor)
plotHistogram_bus(pf_samples[:qd], "qd", "./plots/non-intrusive_2unc"; fignum = 2 + 10, color = buscolor)
plotHistogram_bus(pf_samples[:pg], "pg", "./plots/non-intrusive_2unc"; fignum = 3 + 10, color = buscolor)
plotHistogram_bus(pf_samples[:qg], "qg", "./plots/non-intrusive_2unc"; fignum = 4 + 10, color = buscolor)
plotHistogram_nodal(pf_samples[:e], "e", "./plots/non-intrusive_2unc"; figbum = 5 + 10, color = nodecolor)
plotHistogram_nodal(pf_samples[:f], "f", "./plots/non-intrusive_2unc"; figbum = 6 + 10, color = nodecolor)


### POST PROCESSING ###
width, height, color = "3.9cm", "2.75cm", "black!40"
files_to_save = Dict(:v => Dict("name" => "voltage_magnitude",
        "color" => color,
        "width" => "2.95cm",
        "height" => height),
    :θ => Dict("name" => "voltage_angle",
        "color" => color,
        "width" => "2.95cm",
        "height" => height),
    :pg => Dict("name" => "active_power",
        "color" => color,
        "width" => width,
        "height" => height),
    :qg => Dict("name" => "reactive_power",
        "color" => color,
        "width" => width,
        "height" => height),
    :pd => Dict("name" => "active_power_demand",
        "color" => color,
        "width" => "4.7cm",
        "height" => height),
    :qd => Dict("name" => "reactive_power_demand",
        "color" => color,
        "width" => "4.7cm",
        "height" => height),
    :i => Dict("name" => "current_magnitude",
        "color" => color,
        "width" => width,
        "height" => height))
#  :pl_t => Dict("name"=>"active_power_to",
# 			"color"=>color,
# 			"width"=>width,
# 			"height"=>height),
#  :ql_t => Dict("name"=>"reactive_power_to",
# 			"color"=>color,
# 			"width"=>width,
# 			"height"=>height)

# Store bus and branch data as csv and tikz files
println()
createCSV(files_to_save, pf_samples)
println("CSV data saved to /csv.")
createTikz(files_to_save, pf_samples, "../csv/")
println("Tikz files saved to /tikz.")