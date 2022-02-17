using PowerFlowUnderUncertainty, LinearAlgebra, JuMP, Ipopt, DelimitedFiles, JLD, Distributions
using Random
Random.seed!(1234)

include("powersystem.jl")

numSamples = 100
numUnc = 2
sys = setupPowerSystem()
μ, σ, w = [2.1, 3.2], [0.3, 0.4], [0.3, 0.7]
@assert sum(w) == 1 "The weights do not sum to one."

### Monte Carlo reference for stochastic power flow with 2 uncertainties ###
## Take samples of power values, compute PF, compute moments.

println("Setting up PF model.")
pf = Model(with_optimizer(Ipopt.Optimizer))
set_optimizer_attribute(pf, "print_level", 2) # set verbosity of Ipopt output. Default is 5.
addCoreDeterministic!(pf, sys)
addPVBusDeterministic!(pf)

### Perform full MC simulation ##
pf_samples = Dict(:pg => Array{Float64}(undef, 2, 0),
    :qg => Array{Float64}(undef, 2, 0),
    :e => Array{Float64}(undef, 4, 0),
    :f => Array{Float64}(undef, 4, 0),
    :i_re => Array{Float64}(undef, 5, 0),
    :i_im => Array{Float64}(undef, 5, 0))

## Model
# Input:  x - vector of sampled values for active power of PQ buses.
# Return: dict of PF outputs (buses and line currents).
function model(x::AbstractVector)
    p1, q1 = x[1], x[1] * 0.85
    p2, q2 = x[2], x[2] * 0.85
    sys[:P][1] = p1
    sys[:Q][1] = q1
    sys[:P][2] = p2
    sys[:Q][2] = q2
    resetPowerFlowConstraint!(pf, sys) # initialize pf model with new sample
    optimize!(pf) # actual model function

    d = Dict(:pg => value.(pf[:pg]),
        :qg => value.(pf[:qg]),
        :e => value.(pf[:e]),
        :f => value.(pf[:f]))

    currents = computeLineCurrentsDeterministic(pf, sys)

    merge!(d, currents)
end

# Generate samples for load buses
samples = [sampleFromGaussianMixture(numSamples, μ, σ, w) for i in 1:numUnc]
X = hcat(samples...)
pf_samples[:pd] = X'
pf_samples[:qd] = X' * 0.85

# Evaluate PF: Y = model.(X)
println("Running $numSamples deterministic PF calculations (model evalutations)...")
i = 1
for x in eachrow(X)
    res = model(x)
    pf_samples[:pg] = hcat(pf_samples[:pg], res[:pg])
    pf_samples[:qg] = hcat(pf_samples[:qg], res[:qg])
    pf_samples[:e] = hcat(pf_samples[:e], res[:e])
    pf_samples[:f] = hcat(pf_samples[:f], res[:f])
    pf_samples[:i_re] = hcat(pf_samples[:i_re], res[:i_re])
    pf_samples[:i_im] = hcat(pf_samples[:i_im], res[:i_im])
    i % 500 == 0 ? println("iteration $i") : nothing
    global i += 1
end

# Additional polar values of current and voltage
computePolarValues!(pf_samples)
println("\n Monter Calro model evaluations for $numSamples samples:")
display(pf_samples)
println()



### Compute and store first three Moments ###
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

# Store moments
f_moms = "coefficients/MC_moments.jld"
save(f_moms, "moments", moments)
println("Monte Carlo moments data saved to $f_coeff.\n")



### Plotting ###
buscolor = "red"
nodecolor = "blue"
plotHistogram_bus(pf_samples[:pd], "pd", "./plots/2unc_Monte_Carlo"; fignum = 1 + 10, color = buscolor)
plotHistogram_bus(pf_samples[:qd], "qd", "./plots/2unc_Monte_Carlo"; fignum = 2 + 10, color = buscolor)
plotHistogram_bus(pf_samples[:pg], "pg", "./plots/2unc_Monte_Carlo"; fignum = 3 + 10, color = buscolor)
plotHistogram_bus(pf_samples[:qg], "qg", "./plots/2unc_Monte_Carlo"; fignum = 4 + 10, color = buscolor)
plotHistogram_nodal(pf_samples[:e], "e", "./plots/2unc_Monte_Carlo"; figbum = 5 + 10, color = nodecolor)
plotHistogram_nodal(pf_samples[:f], "f", "./plots/2unc_Monte_Carlo"; figbum = 6 + 10, color = nodecolor)



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