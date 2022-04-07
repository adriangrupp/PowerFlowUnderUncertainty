using PowerFlowUnderUncertainty, PowerModels, LinearAlgebra, Ipopt, JuMP, JLD

"""
30 Bus net: Monte Carlo reference for stochastic power flow
1 Uncertainty: Bus 8 (Load 5)
Take samples of power values, compute PF, perform regression for all needed variables.
"""

caseFile = "case30.m"
numSamples = 10000
Nunc = 1
maxDeg = 4
postProcessing = true

println("\n\t\t===== Stochastic Power Flow: 30 Bus case, 1 Uncertainty, Monte Carlo Simulation =====\n")

## Read case file, initialize network uncertainties and corresponding values
include("init_ni.jl")
network_data = readCaseFile(caseFile)
sys = parseNetworkData(network_data)
# Define uncertain buses
p = sys[:Pd][5]
q = sys[:Qd][5]
unc = initUncertainty_Nu(p, q, numSamples)

# Dict for results on bus and branch parameters
pf_samples = Dict(:pg => Array{Float64}(undef, sys[:Ng], 0),
    :qg => Array{Float64}(undef, sys[:Ng], 0),
    :e => Array{Float64}(undef, sys[:Nbus], 0),
    :f => Array{Float64}(undef, sys[:Nbus], 0),
    :i_re => Array{Float64}(undef, sys[:Nline], 0),
    :i_im => Array{Float64}(undef, sys[:Nline], 0)
)

## Initialize solver
solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 2)

## Retreive samples for initialized uncertainties
pf_samples[:pd] = unc[:samples_bus][:, 1]' # Transpose since we store as rows
pf_samples[:qd] = unc[:samples_bus][:, 2]'

## Execute the model for all samples
println("Running $numSamples deterministic PF calculations (model evalutations)...")
i = 1
for x in eachrow(unc[:samples_bus])
    i % 1000 == 0 ? println("Iteration $i") : nothing
    global i += 1

    network_data["load"]["5"]["pd"] = x[1] # bus 8 / load 5
    network_data["load"]["5"]["qd"] = x[2] # bus 8 / load 5

    res = runPfModel(network_data, solver)
    currents = computeLineCurrentsDeterministic(res[:e], res[:f], sys)
    merge!(res, currents)

    pf_samples[:pg] = hcat(pf_samples[:pg], res[:pg])
    pf_samples[:qg] = hcat(pf_samples[:qg], res[:qg])
    pf_samples[:e] = hcat(pf_samples[:e], res[:e])
    pf_samples[:f] = hcat(pf_samples[:f], res[:f])
    pf_samples[:i_re] = hcat(pf_samples[:i_re], res[:i_re])
    pf_samples[:i_im] = hcat(pf_samples[:i_im], res[:i_im])
end

## Additional polar values of current and voltage
computePolarValues!(pf_samples)

println("\n Monter Calro model evaluations for $numSamples samples:")
display(pf_samples)
println()

### Compute and store first two moments ###
moments = Dict{Symbol,Matrix{Float64}}()
for (key, samples) in pf_samples
    let moms = Array{Float64}(undef, 0, 2)
        for row in eachrow(samples)
            mean_mc, std_mc = mean(row), std(row) #, skewness(row)
            moms = vcat(moms, [mean_mc std_mc])
        end
        moments[key] = moms
    end
end

println("\n Moments of bus parameters:")
display(moments)
println()

## Store moments
f_moms = "coefficients/SPF_MC_moments.jld"
save(f_moms, "moments", moments)
println("Monte Carlo moments data saved to $f_moms.\n")


### POST PROCESSING ###
if postProcessing
    mycolor = "red"
    plotHistogram_6in9(pf_samples[:pg], "pg", "./plots/monte-carlo"; fignum=1 + 10, color=mycolor)
    plotHistogram_6in9(pf_samples[:qg], "qg", "./plots/monte-carlo"; fignum=2 + 10, color=mycolor)

    plotHistogram_9in9(pf_samples[:e][1:9, :], "e1", "./plots/monte-carlo"; fignum=3 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:e][10:18, :], "e2", "./plots/monte-carlo"; fignum=4 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:e][19:27, :], "e3", "./plots/monte-carlo"; fignum=5 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:e][28:30, :], "e4", "./plots/monte-carlo"; fignum=6 + 10, color=mycolor)

    plotHistogram_9in9(pf_samples[:f][1:9, :], "f1", "./plots/monte-carlo"; fignum=7 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:f][10:18, :], "f2", "./plots/monte-carlo"; fignum=8 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:f][19:27, :], "f3", "./plots/monte-carlo"; fignum=9 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:f][28:30, :], "f4", "./plots/monte-carlo"; fignum=10 + 10, color=mycolor)

    plotHistogram_unc(pf_samples[:pd][:], pf_samples[:qd][:], ["pd", "qd"], "./plots/monte-carlo"; fignum=0 + 10, color=mycolor)
end