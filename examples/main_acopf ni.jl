using PowerFlowUnderUncertainty, LinearAlgebra, JuMP, Ipopt, DelimitedFiles
include("powersystem.jl")
include("init_sparse.jl")
# close("all")

### Non-intrusive PCE ###
## Take samples of power values, compute PF, perform regression for all needed variables.
numSamples = 30
maxDegree = deg

# unconstrained problem
println("Setting up OPF model.")
opf = Model(with_optimizer(Ipopt.Optimizer, max_iter = 1000))
addCoreDeterministic!(opf, sys)
addConsistencyConstraints!(opf) # TODO
addCost!(opf, sys, unc) # TODO

# Wrapper function for NI-algo. For each sample of the PQ bus values evaluate the OPF once
# Input:  x - sampled value for active power of PQ bus.
# Return: dict of PF outputs.
function model(x)
    p, q = x, x * 0.85
    sys[:P][1] = p
    sys[:Q][1] = q
    resetPowerFlowConstraint!(pf, sys)
    optimize!(pf) # model function

    return Dict(:pg => value.(pf[:pg]),
        :qg => value.(pf[:qg]),
        :e => value.(pf[:e]),
        :f => value.(pf[:f]))
end

# Create samples
X = sampleFromGaussianMixture(numSamples, μ, σ, w)

# Perform full non-intrusive computations
busRes = Dict(:pg => Array{Float64}(undef, 2, 0), # TODO?
    :qg => Array{Float64}(undef, 2, 0),   # TODO?
    :e => Array{Float64}(undef, 4, 0),
    :f => Array{Float64}(undef, 4, 0))

# Y = model.(X)
for x in X
    res = model(x)
    busRes[:pg] = hcat(busRes[:pg], res[:pg])
    busRes[:qg] = hcat(busRes[:qg], res[:qg])
    busRes[:e] = hcat(busRes[:e], res[:e])
    busRes[:f] = hcat(busRes[:f], res[:f])
end

# Perform the actual regression for PCE coefficients on pd, qd, e and f
pce = computeNonIntrusiveCoefficients(X, busRes, maxDegree, unc) #TODO? for new variables

# Get PCE of currents, branch flows and demands
opf_state = getGridState(opf, sys, unc) #TODO
# Sample for parameter values, using their PCE representations
opf_samples = generateSamples(ξ, opf_state, sys, unc)
display(opf_samples) # Debugging

# constrained problem
# opf_con = Model(with_optimizer(Ipopt.Optimizer, max_iter = 1000))
# addCore!(opf_con, sys, unc)
# addConsistencyConstraints!(opf_con)
# addGenerationConstraint!(2, :pg, opf_con, sys, unc)
# addGenerationConstraint!(1, :qg, opf_con, sys, unc)
# addCost!(opf_con, sys, unc)
# optimize!(opf_con)

# opf_con_state = getGridState(opf_con, sys, unc)
# opf_con_samples = generateSamples(ξ, opf_con_state, sys, unc)

########################################################################
##### POST PROCESSING #####
########################################################################

mycolor = "red"
plotHistogram_gen(opf_samples[:pg], "pg"; fignum = 1, color = mycolor)
plotHistogram_gen(opf_samples[:qg], "qg"; fignum = 2, color = mycolor)
plotHistogram_nodal(opf_samples[:v], "v"; fignum = 3, color = mycolor)
plotHistogram_nodal(opf_samples[:θ], "θ"; fignum = 4, color = mycolor)
plotHistogram_branch(opf_samples[:i], "i"; fignum = 5, color = mycolor)
plotHistogram_branch(opf_samples[:pl_t], "pl_t"; fignum = 6, color = mycolor)

width, height, color, compcolor = "3.9cm", "2.75cm", "blue!40", "gray!20"
files_to_save = Dict(:v => Dict("name" => "voltage_magnitude",
        "color" => color,
        "compcolor" => compcolor,
        "width" => "2.95cm",
        "height" => height),
    :θ => Dict("name" => "voltage_angle",
        "color" => color,
        "compcolor" => compcolor,
        "width" => "2.95cm",
        "height" => height),
    :pg => Dict("name" => "active_power",
        "color" => color,
        "compcolor" => compcolor,
        "width" => "2.95cm",
        "height" => height),
    :qg => Dict("name" => "reactive_power",
        "color" => color,
        "compcolor" => compcolor,
        "width" => "2.95cm",
        "height" => height),
    :pd => Dict("name" => "active_power_demand",
        "color" => color,
        "compcolor" => compcolor,
        "width" => "4.7cm",
        "height" => height),
    :qd => Dict("name" => "reactive_power_demand",
        "color" => color,
        "compcolor" => compcolor,
        "width" => "4.7cm",
        "height" => height),
    :i => Dict("name" => "current_magnitude",
        "color" => color,
        "compcolor" => compcolor,
        "width" => width,
        "height" => height),
    :pl_t => Dict("name" => "active_power_to",
        "color" => color,
        "compcolor" => compcolor,
        "width" => width,
        "height" => height),
    :ql_t => Dict("name" => "reactive_power_to",
        "color" => color,
        "compcolor" => compcolor,
        "width" => width,
        "height" => height))

createCSV(files_to_save, opf_samples)
createTikz(files_to_save, opf_samples, "", "../ppf/")

# constrained problem
##################################################
# mycolor = "green"
# plotHistogram_gen(opf_con_samples[:pg], "pg"; fignum = 1, color = mycolor)
# plotHistogram_gen(opf_con_samples[:qg], "qg"; fignum = 2, color = mycolor)
# plotHistogram_nodal(opf_con_samples[:v], "v"; fignum = 3, color = mycolor)
# plotHistogram_nodal(opf_con_samples[:θ], "θ"; fignum = 4, color = mycolor)
# plotHistogram_branch(opf_con_samples[:i], "i"; fignum = 5, color = mycolor)
# plotHistogram_branch(opf_con_samples[:pl_t], "pl_t"; fignum = 6, color = mycolor)

# width, height, color = "3.9cm", "2.75cm", "red!40"
# files_to_save = Dict(:v => Dict("name"=>"voltage_magnitude",
# 								"color"=>color,
# 								"compcolor"=>"blue!20",
# 								"width"=>"2.95cm",
# 								"height"=>height),
# 					 :θ => Dict("name"=>"voltage_angle",
# 								"color"=>color,
# 								"compcolor"=>"blue!20",
# 								"width"=>"2.95cm",
# 								"height"=>height),
# 					 :pg => Dict("name"=>"active_power",
# 								"color"=>color,
# 								"compcolor"=>"blue!20",
# 								"width"=>"2.95cm",
# 								"height"=>height),
# 					 :qg => Dict("name"=>"reactive_power",
# 								"color"=>color,
# 								"compcolor"=>"blue!20",
# 								"width"=>"2.95cm",
# 								"height"=>height),
# 					 :pd => Dict("name"=>"active_power_demand",
# 								"color"=>color,
# 								"compcolor"=>"blue!20",
# 								"width"=>"4.7cm",
# 								"height"=>height),
# 					 :qd => Dict("name"=>"reactive_power_demand",
# 								"color"=>color,
# 								"compcolor"=>"blue!20",
# 								"width"=>"4.7cm",
# 								"height"=>height),
# 					 :i => Dict("name"=>"current_magnitude",
# 								"color"=>color,
# 								"compcolor"=>"blue!20",
# 								"width"=>width,
# 								"height"=>height),
# 					 :pl_t => Dict("name"=>"active_power_to",
# 								"color"=>color,
# 								"compcolor"=>"blue!20",
# 								"width"=>width,
# 								"height"=>height),
# 					 :ql_t => Dict("name"=>"reactive_power_to",
# 								"color"=>color,
# 								"compcolor"=>"blue!20",
# 								"width"=>width,
# 								"height"=>height))

# createCSV(files_to_save,opf_con_samples)
# createTikz(files_to_save,opf_con_samples,"","../acopf/")
