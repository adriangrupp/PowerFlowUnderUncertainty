using PowerFlowUnderUncertainty, PowerModels, LinearAlgebra, Ipopt, JuMP, TimerOutputs, JLD2

"""
30 Bus net: Monte Carlo reference for stochastic power flow
10 Uncertainties: Bus 3(l2), 4(l3), 7(l4), 8(l5), 10(l6), 17(l11), 19(l13), 29(l19), 30(l20)  and Bus 13(g6).
Take samples of power values, compute PF, perform regression for all needed variables.
"""

caseFile = "case30.m"
numSamples = 10000
nUnc = 10
maxDeg = 3
postProcessing = false

println("\n\t\t===== Stochastic Power Flow: 30 Bus case, 10 Uncertainties, Monte Carlo Simulation =====\n")

## Read case file, initialize network uncertainties and corresponding values
include("init_ni.jl")
network_data = readCaseFile(caseFile)
sys = parseNetworkData(network_data)
# Define uncertain buses
p = [sys[:Pd][2],
    sys[:Pd][3],
    sys[:Pd][4],
    sys[:Pd][5],
    sys[:Pd][6],
    sys[:Pd][11],
    sys[:Pd][13],
    sys[:Pd][19],
    sys[:Pd][20],
    sys[:Pg][6]]
q = [sys[:Qd][2],
    sys[:Qd][3],
    sys[:Qd][4],
    sys[:Qd][5],
    sys[:Qd][6],
    sys[:Qd][11],
    sys[:Qd][13],
    sys[:Qd][19],
    sys[:Qd][20],
    sys[:Qg][6]]
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
pf_samples[:pd] = unc[:samples_bus][:, 1:nUnc]' # Transpose since we store as rows
pf_samples[:qd] = unc[:samples_bus][:, nUnc+1:end]'

## Timer for profiling
to = TimerOutput()

## Execute the model for all samples
println("Running $numSamples deterministic PF calculations (model evalutations)...")
i = 1
@timeit to "Model evaluations" begin
    for x in eachrow(unc[:samples_bus]) # each row is a sample set
        i % 1000 == 0 ? println("Iteration $i") : nothing
        global i += 1

        network_data["load"]["2"]["pd"] = x[1]      # bus 3 / load 2
        network_data["load"]["2"]["qd"] = x[nUnc+1] # second half of matrix are q values
        network_data["load"]["3"]["pd"] = x[2]
        network_data["load"]["3"]["qd"] = x[nUnc+2]
        network_data["load"]["4"]["pd"] = x[3]
        network_data["load"]["4"]["qd"] = x[nUnc+3]
        network_data["load"]["5"]["pd"] = x[4]      # bus 8 / load 5
        network_data["load"]["5"]["qd"] = x[nUnc+4]
        network_data["load"]["6"]["pd"] = x[5]
        network_data["load"]["6"]["qd"] = x[nUnc+5]
        network_data["load"]["11"]["pd"] = x[6]
        network_data["load"]["11"]["qd"] = x[nUnc+6]
        network_data["load"]["13"]["pd"] = x[7]
        network_data["load"]["13"]["qd"] = x[nUnc+7]
        network_data["load"]["19"]["pd"] = x[8]
        network_data["load"]["19"]["qd"] = x[nUnc+8]
        network_data["load"]["20"]["pd"] = x[9]
        network_data["load"]["20"]["qd"] = x[nUnc+9]
        network_data["gen"]["6"]["pg"] = x[10]      # bus 13 / gen 6, no q values vor generators

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
end
println("Finished.")

## Additional polar values of current and voltage
computePolarValues!(pf_samples)

println("\n Monter Calro model evaluations for $numSamples samples:")
display(pf_samples)
println()

### Compute and store first two moments ###
@timeit to "Moments calculation" begin
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
end

println("\n Moments of bus parameters:")
display(moments)
println()

## Store moments
f_moms = "coefficients/SPF_10u_MC_moments.jld2"
save(f_moms, "moments", moments)
println("Monte Carlo moments data saved to $f_moms.\n")

## Show timing stats
println("Timing resutls:")
show(to)


### POST PROCESSING ###
if postProcessing
    mycolor = "red"
    plotHistogram_6in9(pf_samples[:pg], "pg", "./plots/10u_monte-carlo"; fignum=1 + 10, color=mycolor)
    plotHistogram_6in9(pf_samples[:qg], "qg", "./plots/10u_monte-carlo"; fignum=2 + 10, color=mycolor)

    plotHistogram_9in9(pf_samples[:e][1:9, :], "e1", "./plots/10u_monte-carlo"; fignum=3 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:e][10:18, :], "e2", "./plots/10u_monte-carlo"; fignum=4 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:e][19:27, :], "e3", "./plots/10u_monte-carlo"; fignum=5 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:e][28:30, :], "e4", "./plots/10u_monte-carlo"; fignum=6 + 10, color=mycolor)

    plotHistogram_9in9(pf_samples[:f][1:9, :], "f1", "./plots/10u_monte-carlo"; fignum=7 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:f][10:18, :], "f2", "./plots/10u_monte-carlo"; fignum=8 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:f][19:27, :], "f3", "./plots/10u_monte-carlo"; fignum=9 + 10, color=mycolor)
    plotHistogram_9in9(pf_samples[:f][28:30, :], "f4", "./plots/10u_monte-carlo"; fignum=10 + 10, color=mycolor)
    # Plot P & Q of uncertainties
    plotHistogram_2unc(pf_samples[:pd], pf_samples[:qd], ["pd_8", "qd_8", "pg_13", "qg_13"], "./plots/10u_monte-carlo"; fignum=0 + 10, color=mycolor)
end