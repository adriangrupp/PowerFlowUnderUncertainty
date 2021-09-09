using PowerFlowUnderUncertainty, LinearAlgebra, JuMP, Ipopt, DelimitedFiles
include("powersystem.jl")
include("init_sparse.jl")

### Non-intrusive PCE ###
## Take samples of power values, compute PF, perform regression for all needed variables.
numSamples = 100
maxDegree = 4

println("Setting up OPF model.")
pf = Model(with_optimizer(Ipopt.Optimizer))
addCoreDeterministic!(pf,sys)

## Model: Wrapper function for NI-algo.
# Input: x - sampled value for active power of PQ bus.
# Returns: dict of PF outputs.
function model(x)
	p, q = x, x * 0.85
	sys[:P][1] = p
	sys[:Q][1] = q
	resetPowerFlowConstraint!(pf, sys)
	optimize!(pf) # model function
	
	return Dict(:pg=>value.(pf[:pg]),
				:qg=>value.(pf[:qg]),
				:e=> value.(pf[:e]),
				:f=> value.(pf[:f]))
end

busRes = Dict(:pg=>Array{Float64}(undef, 2, 0),
			  :qg=>Array{Float64}(undef, 2, 0),
			  :e=> Array{Float64}(undef, 4, 0),
			  :f=> Array{Float64}(undef, 4, 0))

# Full non-intrusive computations
X = sampleFromGaussianMixture(numSamples,μ,σ,w)
# Y = model.(X)
for x in X
	res = model(x)
	busRes[:pg] = hcat(busRes[:pg], res[:pg])
	busRes[:qg] = hcat(busRes[:qg], res[:qg])
	busRes[:e]  = hcat(busRes[:e], res[:e])
	busRes[:f]  = hcat(busRes[:f], res[:f])
end

pce = computeNonIntrusiveCoefficients(X, busRes, maxDegree, unc)

pf_state = getGridStateNonintrusive(pce, pf, sys, unc)
display(pf_state)

pf_samples = generateSamples(ξ,pf_state,sys,unc) # sample for results, using the built PCE model
display(pf_samples)


### Sparse PCE ###
# Solve sparse PCE with Power flow as model
# sampleFun(sampleSize) = sampleFromGaussianMixture(sampleSize,μ,σ,w)
# pce, Ap, p, R², Q² = sparsePCE(unc[:opq], model, sampleFun)


### Plotting ###
mycolor = "red"
plotHistogram_gen(pf_samples[:pd], "pd"; fignum = 1+10, color = mycolor)
plotHistogram_gen(pf_samples[:qd], "qd"; fignum = 2+10, color = mycolor)
plotHistogram_gen(pf_samples[:pg], "pg"; fignum = 3+10, color = mycolor)
plotHistogram_gen(pf_samples[:qg], "qg"; fignum = 4+10, color = mycolor)
plotHistogram_nodal(pf_samples[:e], "e"; figbum = 5+10, color = mycolor)
plotHistogram_nodal(pf_samples[:f], "f"; figbum = 6+10, color = mycolor)


### POST PROCESSING ###
width, height, color = "3.9cm", "2.75cm", "black!40"
files_to_save = Dict(:v => Dict("name"=>"voltage_magnitude",
								"color"=>color,
								"width"=>"2.95cm",
								"height"=>height),
					 :θ => Dict("name"=>"voltage_angle",
								"color"=>color,
								"width"=>"2.95cm",
								"height"=>height),
					 :pg => Dict("name"=>"active_power",
								"color"=>color,
								"width"=>width,
								"height"=>height),
					 :qg => Dict("name"=>"reactive_power",
								"color"=>color,
								"width"=>width,
								"height"=>height),
					 :pd => Dict("name"=>"active_power_demand",
								"color"=>color,
								"width"=>"4.7cm",
								"height"=>height),
					 :qd => Dict("name"=>"reactive_power_demand",
								"color"=>color,
								"width"=>"4.7cm",
								"height"=>height),
					 :i => Dict("name"=>"current_magnitude",
								"color"=>color,
								"width"=>width,
								"height"=>height))
					#  :pl_t => Dict("name"=>"active_power_to",
					# 			"color"=>color,
					# 			"width"=>width,
					# 			"height"=>height),
					#  :ql_t => Dict("name"=>"reactive_power_to",
					# 			"color"=>color,
					# 			"width"=>width,
					# 			"height"=>height)

createCSV(files_to_save,pf_samples)
createTikz(files_to_save,pf_samples,"")



