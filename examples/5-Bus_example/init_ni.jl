# Initialize for (sparse) non-intrusive PCE - global parameters, uncertainties, polynomial basis, and samples
using Random, PolyChaos
Random.seed!(1234)

"""
Read in network data
"""
function readCaseFlie(caseFile)
    PowerModels.logger_config!("error") # Suppress warnings on phase angles
    network_data = PowerModels.parse_file(caseFile)
    network_data = make_basic_network(network_data) # no dc lines, switches, inactive components, ...
    # print_summary(network_data)
end

"""
Parse network data
"""
function parseNetworkData(network_data)
    N_bus = length(network_data["bus"])
    N_load = length(network_data["load"])
    N_gen = length(network_data["gen"])
    N_branch = length(network_data["branch"])

    incidence = calc_basic_incidence_matrix(network_data) * -1 # we use inverse notation for directions
    Zbr = calc_basic_branch_series_impedance(network_data)
    Ybr = Diagonal(inv.(Zbr)) # branch admittance

    # P and Q values for load buses
    # load_idx = getLoadIndices(network_data)
    Pload, Qload = getLoadPQ(network_data)

    # P and Q values for bus injection of all buses (these are not nec. PQ buses or contain loads)
    # injection = calc_basic_bus_injection(network_data)
    # P = real(injection)
    # Q = imag(injection)

    # generator cost
    costquad, costlin = getCost(network_data)

    sys = Dict(
        :Nbus => N_bus,
        :Nd => N_load,
        :Ng => N_gen,
        :Nline => N_branch,
        :A => incidence,
        :Ybr => Ybr,
        :Pd => Pload,
        :Qd => Qload,
        :costquad => costquad,
        :costlin => costlin
    )
end

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

## Managing of the OPF model. Wrapper function for NI-algo. Currently hard coded for bus 8 / load 5
# Input:  x - sampled value for active power of PQ bus.
# Return: dict of OPF outputs.
function runOpfModel(network_data, solver)
    pf = PowerModels.run_opf(network_data, ACRPowerModel, solver)

    # check if solution was feasible
    status = pf["termination_status"]
    status == OPTIMAL || status == LOCALLY_SOLVED ? nothing : error("Potentially no solution found: ", status)

    pg, qg = getGenPQResult(pf)
    vr, vi = getVoltageResult(pf)

    return Dict(:pg => pg,
        :qg => qg,
        :e => vr, # we use e and f as convention for real and imaginary voltage parts
        :f => vi)
end