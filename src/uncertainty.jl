export  sampleFromGaussianMixture,
        setupUncertainty,
        setupUncertaintySparse,
        computeNonIntrusiveCoefficients,
        generateSamples

function ρ_gauss(x,μ,σ)
    1 / sqrt(2*π*σ^2) * exp(-(x - μ)^2 / (2σ^2))
end

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

function setupUncertaintySparse(μ::Vector,σ::Vector,w::Vector,n::Int,deg::Int)
    @assert length(μ) == length(σ) == length(w) "inconsistent lengths of μ and σ"
    ρ(x) = sum( w[i]*ρ_gauss(x,μ[i],σ[i]) for i in 1:length(w) )
    meas = Measure("my_GaussMixture", ρ, (-Inf,Inf), false, Dict(:μ=>μ,:σ=>σ,:w=>w)) # build measure
    opq = OrthoPoly("my_op",deg,meas;Nquad=150,Nrec = 5*deg, discretization=stieltjes) # construct orthogonal polynomial
    showbasis(opq,digits=2) # in case you wondered
    # PCE of demands
    pd = zeros(n,deg+1)
    pd[1, [1,2]] = calculateAffinePCE(opq) # random load (bus 2)
    pd[2, 1] = 1.2 # deterministic load (bus 4)
    qd = 0.85 * copy(pd)

    return Dict(:opq=>opq,
                :dim=>opq.deg+1,
                :pd=>pd,
                :qd=>qd)
end


# Compute the non-intrusive pce coefficients by least squares regression
function computeNonIntrusiveCoefficients(X::Vector, busRes::Dict, maxDeg::Int, unc::Dict)
    # Evaluate polynomial basis up to maxDegree for all X samples
    Φ = [ evaluate(j, X[i], unc[:opq]) for i = 1:length(X), j = 0:maxDeg ]

    dim = unc[:dim]
    pg = Array{Float64}(undef, 0, dim)
	qg = Array{Float64}(undef, 0, dim)
	e =  Array{Float64}(undef, 0, dim)
	f =  Array{Float64}(undef, 0, dim)

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

    return Dict(:pg=>pg, :qg=>qg, :e=> e, :f=> f)
end


function sampleFromGaussianMixture(n::Int,μ::Vector{},σ::Vector{},w::Vector{})
    X = Float64[]
    for i in 1:n
        k = findfirst(x -> x > rand(), cumsum(w))
        push!(X, μ[k] + σ[k]*randn())
    end
    return X
end

function generateSamples(x,d_in::Dict,sys::Dict,unc::Dict)
    Φ = evaluate(x,unc[:opq])
    d_out = Dict{Symbol,Matrix{Float64}}()
    for (key, value) in d_in
        d_out[key] = (Φ*value')'
    end
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
    X = Matrix{Float64}(undef,size(x,1),1)
    X[:] = [ sum( p[i]^2*sys[:costquad][i] + p[i]*sys[:costlin][i] for i in 1:sys[:Ng] ) for p in eachcol(d_out[:pg]) ]
    d_out[:cost] = X
    return d_out
end
