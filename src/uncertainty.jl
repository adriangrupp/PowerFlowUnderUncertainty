export  sampleFromGaussianMixture,
        setupUncertainty,
        setupUncertaintySparse,
        setupUncertaintyMulti,
        computeCoefficientsNI,
        computeCoefficientsSparse,
        generateSamples,
        computePolarValues!,
        computeMoments

function ρ_gauss(x,μ,σ)
    1 / sqrt(2*π*σ^2) * exp(-(x - μ)^2 / (2σ^2))
end

# Setup Intrusive univariate uncertainty
function setupUncertainty(μ::Vector,σ::Vector,w::Vector,n::Int,deg::Int)
    @assert length(μ) == length(σ) == length(w) "inconsistent lengths of μ and σ"
    ρ(x) = sum( w[i]*ρ_gauss(x,μ[i],σ[i]) for i in 1:length(w) )
    meas = Measure("my_GaussMixture", ρ, (-Inf,Inf), false, Dict(:μ=>μ,:σ=>σ,:w=>w)) # build measure
    opq = OrthoPoly("my_op",deg,meas;Nquad=150,Nrec = 5*deg, discretization=stieltjes) # construct orthogonal polynomial
    showbasis(opq,digits=2) # in case you wondered

    pd = zeros(n,deg+1) # demands
    pd[1, [1,2]] = calculateAffinePCE(opq) # random load (bus 2)
    pd[2, 1] = 1.2 # deterministic load (bus 4)
    qd = 0.85 * copy(pd)

    return Dict(:opq=>opq,
                :T2=>Tensor(2,opq),
                :T3=>Tensor(3,opq),
                :T4=>Tensor(4,opq),
                :pd=>pd,
                :qd=>qd,
                :dim=>size(pd,2))
end

# Setup Non-intrusive univariate Gaussian mixture uncertainty
function setupUncertaintySparse(μ::Vector, σ::Vector, w::Vector, n::Int, deg::Int)
    @assert length(μ) == length(σ) == length(w) "inconsistent lengths of μ and σ"
    ρ(x) = sum(w[i] * ρ_gauss(x, μ[i], σ[i]) for i in 1:length(w))
    meas = Measure("my_GaussMixture", ρ, (-Inf, Inf), false, Dict(:μ => μ, :σ => σ, :w => w)) # build measure
    opq = OrthoPoly("my_op", deg, meas; Nquad = 150, Nrec = 5 * deg, discretization = stieltjes) # construct orthogonal polynomial
    println("Polynomial basis:")
    showbasis(opq, digits = 2) # in case you wondered
    println()
    # PCE of demands
    pd = zeros(n, deg + 1)
    pd[1, [1, 2]] = calculateAffinePCE(opq) # random load (bus 2)
    pd[2, 1] = 1.2 # deterministic load (bus 4)
    qd = 0.85 * copy(pd)

    return Dict(:opq => opq,
        :dim => opq.deg + 1,
        :pd => pd,
        :qd => qd)
end

# Setup for univariate general uncertainty
function setupUncertainty(deg::Int, op::AbstractOrthoPoly)
    println("Polynomial basis:")
    showbasis(op, digits = 2) # in case you wondered
    println()
    # PCE of demand. Compute affine coefficients for each univariate unvertainty dimension and combine them
    pd = zeros(n, deg + 1)
    pd[1, [1,2]] = calculateAffinePCE(op) # random load 
    qd = copy(pd)

    return Dict(:opq => op,
        :dim => op.deg + 1,
        :pd => pd,
        :qd => qd)
end


# Setup for 2 uncertainties following a Gaussian mixture
function setupUncertaintyMulti(μ::Vector, σ::Vector, w::Vector, deg::Int, numUnc::Int)
    @assert length(μ) == length(σ) == length(w) "inconsistent lengths of μ and σ"
    ρ(x) = sum(w[i] * ρ_gauss(x, μ[i], σ[i]) for i in 1:length(w))
    meas = Measure("GaussMixture", ρ, (-Inf, Inf), false, Dict(:μ => μ, :σ => σ, :w => w)) # build measure
    opq = OrthoPoly("GaussMix_op", deg, meas; Nquad = 150, Nrec = 5 * deg, discretization = stieltjes) # construct orthogonal polynomial
    mop = MultiOrthoPoly([opq for i in 1:numUnc], deg)
    println("Polynomial basis:")
    show(mop) # in case you wondered
    println()
    # PCE of demands. Compute affine coefficients for each univariate unvertainty dimension and combine them
    pd1 = assign2multi(calculateAffinePCE(mop.uni[1]), 1, mop.ind) # random load (bus 2)
    pd2 = assign2multi(calculateAffinePCE(mop.uni[2]), 2, mop.ind) # random load (bus 4)
    pd = vcat(pd1', pd2')
    qd = 0.85 * copy(pd)

    return Dict(:opq => mop,
        :dim => mop.dim,
        :pd => pd,
        :qd => qd)
end

# Setup for multiple uncertainties - for now hard coded for 3 load buses #TODO: fix for pv and finish
function setupUncertaintyMulti(p::Vector, q::Vector, mop::MultiOrthoPoly, deg::Int, numUnc::Int)
    # Sample uncertainties for NI model (samples) and postprocessing (ξ)
    samples = [sampleMeasure(numSamples, op) for op in ops]
    samples = hcat(samples...) # numSamples x nUnc matrix
    ξ = [sampleMeasure(5000, op) for op in ops]  # Evaluation samples for each uncertainty
    ξ = hcat(ξ...) # 5000 x nUnc

    # PCE of demand. Compute affine coefficients for univariate uncertainty. Transform support of distributions
    pd = zeros(nUnc, mop.dim) # matrix storing all pd PCE coefficients
    qd = zeros(nUnc, mop.dim)
    samples_p = zeros(numSamples, nUnc) # matrix storing all p-samples
    samples_q = zeros(numSamples, nUnc)
    δ = 0.15 # deviation from nominal value
    # Iterate for all uncertainties
    for (i, op) in enumerate(ops)
        lp, up = [1 - δ, 1 + δ] * p[i]
        lq, uq = [1 - δ, 1 + δ] * q[i]
        pce_p = supp2pce(α[i], β[i]) * [lp; up]
        pce_q = supp2pce(α[i], β[i]) * [lq; uq]
        pdi = assign2multi(pce_p, i, mop.ind) # Get sparse array assigned with first two pce coefficients to multivariate polynomial
        qdi = assign2multi(pce_q, i, mop.ind)
        pd[i, :] = pdi' # Collect index/coefficient arrays for all uncertainties
        qd[i, :] = qdi'

        samples_p[:, i] = samples[:, i] * (up - lp) .+ lp # Generate samples for each uncertainty
        samples_q[:, i] = samples[:, i] * (uq - lq) .+ lq # Generate samples for each uncertainty
    end

    return Dict(:opq => mop,
        :dim => mop.dim,
        :pd => pd, # PCE of initially uncertain buses TODO: currently this encompasses gens and loads
        :qd => qd,
        :samples_unc => samples,                    # samples for exogenous uncertainties
        :samples_bus => hcat(samples_p, samples_q), # samples for bus uncertainties, first all p then all q values
        :ξ => ξ)
end

"""
Core non-intrusive PCE method. Computes PCE coefficients component-wise with provided method.
Can do either full or sparse non-intrusive PCE.
"""
function computeCoefficients(PCEmethod::Function, X::VecOrMat, busRes::Dict, unc::Dict)
    dim = unc[:dim]
    pceRes = Dict() # PCE coefficients
    pceErr = Dict() # LOO-error of all coefficients

    # Perform least squares regression for all relevant bus parameters and get their pce coefficients
    for (key, res) in busRes
        pce = Array{Float64}(undef, 0, dim) # One matrix per parameter
        err = zeros(0)
        # Iterate for each bus
        for row in eachrow(res)
            coeffs = PCEmethod(unc[:Φ], row)
            pce = vcat(pce, coeffs')
            append!(err, looError(row, unc[:Φ], coeffs))
        end
        pceRes[key] = pce
        pceErr[key] = err
    end

    return pceRes, pceErr
end

"""
Compute PCE coefficients fr a full basis via least squares regression
"""
function computeCoefficientsNI(X::VecOrMat, busRes::Dict, unc::Dict)
    method(Φ, Y) = leastSquares(Φ, Y)
    return computeCoefficients(method, X, busRes, unc)
end

"""
Compute PCE coefficients fr a full basis via least squares regression
"""
function computeCoefficientsSparse(X::VecOrMat, busRes::Dict, unc::Dict; K::Int=2)
    method(Φ, Y) = subspacePursuit(Φ, Y, K)[1] # First return value are the PCE coefficients
    return computeCoefficients(method, X, busRes, unc)
end

function sampleFromGaussianMixture(n::Int,μ::Vector{},σ::Vector{},w::Vector{})
    X = Float64[]
    for i in 1:n
        k = findfirst(x -> x > rand(), cumsum(w))
        push!(X, μ[k] + σ[k]*randn())
    end
    return X
end

# Generate samples from the PCE of all the system parameters.
function generateSamples(x,d_in::Dict,sys::Dict,unc::Dict)
    Φ = evaluate(x,unc[:opq])
    # transpose regression matrix for multivariate bases, because PolyChaos somehow swaps dimensions
    typeof(unc[:opq]) <: MultiOrthoPoly ? Φ = Φ' : nothing 
    d_out = Dict{Symbol,Matrix{Float64}}()
    for (key, value) in d_in
        d_out[key] = (Φ*value')'
    end

    computePolarValues!(d_out)

    X = Matrix{Float64}(undef,size(x,1),1)
    X[:] = [ sum( p[i]^2*sys[:costquad][i] + p[i]*sys[:costlin][i] for i in 1:sys[:Ng] ) for p in eachcol(d_out[:pg]) ]
    d_out[:cost] = X
    return d_out
end

function computePolarValues!(d_out::Dict)
    if haskey(d_out,:i_re) && haskey(d_out,:i_im)
        i = d_out[:i_re] + im * d_out[:i_im]
        d_out[:i] = abs.(i)
        d_out[:i_angle] = angle.(i)
    end
    # s = d_out[:pl_t] + im*d_out[:ql_t]
    # d_out[:s_mag] = abs.(s)
    # d_out[:s_ang] = angle.(s)
    if haskey(d_out,:e) && haskey(d_out,:f) 
        v = d_out[:e] + im * d_out[:f]
        d_out[:v] = abs.(v)
        d_out[:θ] = angle.(v)
    end
end

# Compute the moments from the PCE coefficients for each parameter and bus
function computeMoments(d_in::Dict, unc::Dict)
    moments = Dict{Symbol,Matrix{Float64}}()
    for (key, val) in d_in
        println("Computing moments for $key")
        let moms = Array{Float64}(undef, 0, 2)
            for row in eachrow(val)
                mean, std = PolyChaos.mean(row, unc[:opq]), PolyChaos.std(row, unc[:opq])
                moms = vcat(moms, [mean std])
            end
            moments[key] = moms
        end
    end
    return moments
end