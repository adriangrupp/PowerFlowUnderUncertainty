using PowerFlowUnderUncertainty, LinearAlgebra, JuMP, Ipopt, DelimitedFiles
include("powersystem.jl")
include("init_sparse.jl")

numSamples = 4
maxDegree = 5
pf = Model(with_optimizer(Ipopt.Optimizer))
println("Setting up OPF model.")
  addCoreDeterministic!(pf,sys)
# addVariablesDeterministic!(pf, sys)
# addInitialConditionDeterministic!(pf, sys)
# addSlackBusDeterministic!(pf)
# addPVBusDeterministic!(pf)
# addPowerFlowDeterministic!(pf, sys)
optimize!(pf)

### Non-intrusive PCE ###
# Take initial sample of power values, compute PF, perform regression for all needed variables

# Model: wrapper function for sparse algo. x - sampled value for active power of PQ bus
function model(x)
	p, q = x, x * 0.85
	sys[:P][1] = p
	sys[:Q][1] = q
	resetPowerFlowConstraint!(pf, sys)
	optimize!(pf) # model function
	
	# dict of outputs or matrix
	result = Dict()
	result["pg"] = value.(pf[:pg])
	result["qg"] = value.(pf[:qg])
	result["e"] = value.(pf[:e])
	result["f"] = value.(pf[:f])
	return result
end

# Full non-intrusive computations
X = sampleFromGaussianMixture(numSamples,μ,σ,w)

pg = Array{Float64}(undef, 2, 0)
qg = Array{Float64}(undef, 2, 0)
e = Array{Float64}(undef, 4, 0)
f = Array{Float64}(undef, 4, 0)

# Y = model.(X)
for x in X
	res = model(x)
	global pg = hcat(pg, res["pg"])
	global qg = hcat(qg, res["qg"])
	global e = hcat(e, res["e"])
	global f = hcat(f, res["f"])
end

# Evaluate polynomial basis up to maxDegree
Φ = [ evaluate(j, X[i], unc[:opq]) for i = 1:numSamples, j = 0:maxDegree ]
println(Φ)

# Perform least squares regression for all relevant bus variables and get their pce coefficients
pg_pce = [leastSquares(Φ, row) for row in eachrow(pg)]
qg_pce = [leastSquares(Φ, row) for row in eachrow(qg)]
e_pce = [leastSquares(Φ, row) for row in eachrow(e)]
f_pce = [leastSquares(Φ, row) for row in eachrow(f)]

println(pg_pce)

pf_state = getGridStateDeterministic(pf, sys)


### Sparse PCE ###
# Solve sparse PCE with Power flow as model
# sampleFun(sampleSize) = sampleFromGaussianMixture(sampleSize,μ,σ,w)
# pce, Ap, p, R², Q² = sparsePCE(unc[:opq], model, sampleFun)

# pf_state = getGridStateDeterministic(pf,sys)
# pf_samples = generateSamples(ξ,pf_state,sys,unc) # sample for results, using the built PCE model



### Plotting ###
mycolor = "red"
plotHistogram_gen(pf_samples[:pd], "pd"; fignum = 1+10, color = mycolor)
plotHistogram_gen(pf_samples[:qd], "qd"; fignum = 2+10, color = mycolor)

########################################################################
##### POST PROCESSING #####
########################################################################

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
								"height"=>height),
					 :pl_t => Dict("name"=>"active_power_to",
								"color"=>color,
								"width"=>width,
								"height"=>height),
					 :ql_t => Dict("name"=>"reactive_power_to",
								"color"=>color,
								"width"=>width,
								"height"=>height))

createCSV(files_to_save,pf_samples)
createTikz(files_to_save,pf_samples,"")



