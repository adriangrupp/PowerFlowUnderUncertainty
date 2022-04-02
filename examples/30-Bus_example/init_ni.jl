# Initialize for (sparse) non-intrusive PCE - global parameters, uncertainties, polynomial basis, and samples
using Random, PolyChaos, LinearAlgebra
Random.seed!(1234)

"""
Transformation matrix from β-distribution to affine PCE coefficients
    `[ x_0; x_1 ] = supp2pce(α,β)*[ x_lb; x_ub ]`
"""
function supp2pce(α, β)
    T = 1 / (α + β) * [β α; -1 1]
end

"""
Initialization for 1 uncertainty (β-distributed)
"""
function initUncertainty_1u(p, q)
    #shape parameters
    α = 2
    β = 4.666
    op = Beta01OrthoPoly(maxDeg, α, β)

    println("Polynomial basis of size $(op.deg+1):")
    showbasis(op, digits=2)
    println()

    # PCE of demand. Compute affine coefficients for univariate uncertainty. Transform support of distributions
    pd = zeros(1, maxDeg + 1)
    qd = zeros(1, maxDeg + 1)
    δ = 0.15 # deviation from nominal value
    lp, up = [1 - δ, 1 + δ] * p
    lq, uq = [1 - δ, 1 + δ] * q
    pd[1, [1, 2]] = convert2affinePCE(lp, up, op)
    qd[1, [1, 2]] = convert2affinePCE(lq, uq, op)

    samples = sampleMeasure(numSamples, op)  # NI model samples
    samples_p = samples * (up - lp) .+ lp
    samples_q = samples * (uq - lq) .+ lq
    ξ = sampleMeasure(5000, op)  # Evaluation samples

    return Dict(:opq => op,
        :dim => op.deg + 1,
        :pd => pd,
        :qd => qd,
        :samples_unc => samples,    # samples for exogenous uncertainty
        :samples_bus => hcat(samples_p, samples_q), # samples for bus uncertainty
        :ξ => ξ)
end

"""
 Initialization for 2 uncertainties (2 β-distributions)
"""
function initUncertainty_2u(p::Vector, q::Vector)
    # Setup multivariate basis
    α = [2, 3]
    β = [4.666, 3]
    op1 = Beta01OrthoPoly(maxDeg, α[1], β[1]; Nrec = maxDeg+2)
    op2 = Beta01OrthoPoly(maxDeg, α[2], β[2]; Nrec = maxDeg+2)
    ops = [op1, op2]

    mop = MultiOrthoPoly(ops, maxDeg)
    println("Polynomial basis:")
    show(mop)
    println()

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
        pce_q = supp2pce(α[i], β[i]) * [lq; uq] # need to do this due to 0/0 support
        pce_p = convert2affinePCE(lp, up, op)
    
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
 Initialization for various uncertainties (β-distributions)
"""
function initUncertainty_Nu(p::Vector, q::Vector)
    n = 10

    # Setup multivariate basis
    α = [2, 3, 3, 4, 2, 1, 5, 3, 2, 3]
    β = [4.666, 3, 2, 2, 3, 1, 4, 4, 3, 3]
    ops = [Beta01OrthoPoly(maxDeg, α[i], β[i]; Nrec = maxDeg+2) for i in 1:n]

    mop = MultiOrthoPoly(ops, maxDeg)
    println("Polynomial basis:")
    show(mop)
    println()

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
        pce_p = convert2affinePCE(lp, up, op)
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