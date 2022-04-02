# Initialize for (sparse) non-intrusive PCE - global parameters, uncertainties, polynomial basis, and samples
using Random, PolyChaos
Random.seed!(1234)

"""
Transformation matrix from β-distribution to affine PCE coefficients
    `[ x_0; x_1 ] = supp2pce(α,β)*[ x_lb; x_ub ]`
"""
function supp2pce(α, β)
    T = 1 / (α + β) * [β α; -1 1]
end

"""
Initialization for 1 uncertainty (bus 8/load 5)
"""
function initUncertainty_1u(p, q)
    #shape parameters
    α = 2
    β = 4.666
    op = Beta01OrthoPoly(maxDeg, α, β)

    println("Polynomial basis of size $(op.deg+1):")
    showbasis(op, digits = 2)
    println()

    # PCE of demand. Compute affine coefficients for univariate uncertainty. Transform support of distributions
    pd = zeros(1, maxDeg + 1)
    qd = zeros(1, maxDeg + 1)
    δ = 0.15 # deviation from nominal value
    lp, up = [1 - δ, 1 + δ] * p
    lq, uq = [1 - δ, 1 + δ] * q
    pd[1, [1, 2]] = supp2pce(α, β) * [lp; up]
    qd[1, [1, 2]] = supp2pce(α, β) * [lq; uq]

    samples = sampleMeasure(numSamples, op)  # NI model samples
    samples_p = samples * (up - lp) .+ lp
    samples_q = samples * (uq - lq) .+ lq
    ξ = sampleMeasure(5000, op)  # Evaluation samples

    unc = Dict(:opq => op,
        :dim => op.deg + 1,
        :pd => pd,
        :qd => qd,
        :samples_unc => samples,    # samples for exogenous uncertainty
        :samples_bus => hcat(samples_p, samples_q), # samples for bus uncertainty
        :ξ => ξ)
end