using PowerFlowUnderUncertainty, LinearAlgebra, JuMP, Ipopt, DelimitedFiles
include("powersystem.jl")
include("init.jl")

pf = Model(with_optimizer(Ipopt.Optimizer))
addCore!(pf, sys, unc)
addPVBus!(pf, sys, unc)
optimize!(pf)

pf_state = getGridState(pf, sys, unc)
pce_intr = Dict(
    :qg => pf_state[:qg],
    :pg => pf_state[:pg],
    :e => pf_state[:e],
    :f => pf_state[:f]
)

pf_samples = generateSamples(ξ, pf_state, sys, unc)


### Store PCE coefficients ###
save("coefficients/SPF_intrusive.jld", "pf_state", pf_state)

mycolor = "red"
plotHistogram_bus(pf_samples[:pd], "pd", "./plots/intrusive"; fignum = 1 + 10, color = mycolor)
plotHistogram_bus(pf_samples[:qd], "qd", "./plots/intrusive"; fignum = 2 + 10, color = mycolor)
plotHistogram_bus(pf_samples[:pg], "pg", "./plots/intrusive"; fignum = 3 + 10, color = mycolor)
plotHistogram_bus(pf_samples[:qg], "qg", "./plots/intrusive"; fignum = 4 + 10, color = mycolor)
plotHistogram_nodal(pf_samples[:e], "e", "./plots/intrusive"; figbum = 5 + 10, color = mycolor)
plotHistogram_nodal(pf_samples[:f], "f", "./plots/intrusive"; figbum = 6 + 10, color = mycolor)

########################################################################
##### POST PROCESSING #####
########################################################################

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
        "height" => height),
    :pl_t => Dict("name" => "active_power_to",
        "color" => color,
        "width" => width,
        "height" => height),
    :ql_t => Dict("name" => "reactive_power_to",
        "color" => color,
        "width" => width,
        "height" => height))

createCSV(files_to_save, pf_samples)
println("CSV data saved to /csv.")
createTikz(files_to_save, pf_samples, "")
println("Tikz files saved to /tikz.")
