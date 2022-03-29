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

# Intrusive univariate uncertainty
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

# Non-intrusive univariate Gaussian mixture uncertainty
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

# Univariate general uncertainty
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


# Setup for 2 uncertainties (Gaussian mixture)
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

# Setup for multiple uncertainties - for now hard coded for 3 load buses
function setupUncertaintyMulti(deg::Int, ops::Vector, numUnc::Int)
    mop = MultiOrthoPoly(ops, deg)
    println("Polynomial basis:")
    show(mop) # in case you wondered
    println()
    # PCE of demands. Compute affine coefficients for each univariate unvertainty dimension and combine them
    pd1 = assign2multi(calculateAffinePCE(mop.uni[1]), 1, mop.ind) # random load 
    pd2 = assign2multi(calculateAffinePCE(mop.uni[2]), 2, mop.ind) # random load 
    pd3 = assign2multi(calculateAffinePCE(mop.uni[3]), 3, mop.ind) # random load 
    pd = vcat(pd1', pd2', pd3')
    qd = 0.85 * copy(pd)

    return Dict(:opq => mop,
        :dim => mop.dim,
        :pd => pd,
        :qd => qd)
end


# TODO: merge with sparse method
# Compute the non-intrusive pce coefficients component-wise by least squares regression
function computeCoefficientsNI(X::VecOrMat, busRes::Dict, unc::Dict)
    # Evaluate polynomial basis for all samples. Multiply dispatched for uni and multivar
    Φ = evaluate(X, unc[:opq])
    # Transpose regression matrix for multivariate bases, because PolyChaos somehow swaps dimensions
    typeof(unc[:opq]) <: MultiOrthoPoly ? Φ = Φ' : nothing 
    
    dim = unc[:dim]
    pg = Array{Float64}(undef, 0, dim)
    qg = Array{Float64}(undef, 0, dim)
    e = Array{Float64}(undef, 0, dim)
    f = Array{Float64}(undef, 0, dim)

    # Perform least squares regression for all relevant bus variables and get their pce coefficients
    for row in eachrow(busRes[:pg])
        pg = vcat(pg, leastSquares(Φ, row)')
    end
    for row in eachrow(busRes[:qg])
        qg = vcat(qg, leastSquares(Φ, row)')
    end
    for row in eachrow(busRes[:e])
        e = vcat(e, leastSquares(Φ, row)')
    end
    for row in eachrow(busRes[:f])
        f = vcat(f, leastSquares(Φ, row)')
    end

    return Dict(:pg => pg, :qg => qg, :e => e, :f => f)
end


# Compute the non-intrusive pce coefficients component-wise by sparse regresseion (subspace pursuit)
function computeCoefficientsSparse(X::VecOrMat, busRes::Dict, unc::Dict; K::Int=2)
    # Evaluate polynomial basis for all X samples. Multiply dispatched for uni and multivar
    Φ = evaluate(X, unc[:opq])
    # transpose regression matrix for multivariate bases, because PolyChaos somehow swaps dimensions
    typeof(unc[:opq]) <: MultiOrthoPoly ? Φ = Φ' : nothing

    dim = unc[:dim]
    pg = Array{Float64}(undef, 0, dim)
    qg = Array{Float64}(undef, 0, dim)
    e = Array{Float64}(undef, 0, dim)
    f = Array{Float64}(undef, 0, dim)

    # Perform sparse subspace pursuit regression for all relevant bus variables and get their PCE coefficients
    for row in eachrow(busRes[:pg])
        pg = vcat(pg, subspacePursuit(Φ, row, K)[1]')
    end
    for row in eachrow(busRes[:qg])
        qg = vcat(qg, subspacePursuit(Φ, row, K)[1]')
    end
    for row in eachrow(busRes[:e])
        e = vcat(e, subspacePursuit(Φ, row, K)[1]')
    end
    for row in eachrow(busRes[:f])
        f = vcat(f, subspacePursuit(Φ, row, K)[1]')
    end

    return Dict(:pg => pg, :qg => qg, :e => e, :f => f)
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