# Initialize for (sparse) non-intrusive PCE - global parameters, uncertainties, polynomial basis, and samples
using Random, PolyChaos

Random.seed!(1234)

# Parse network data
PowerModels.logger_config!("error")
network_data = PowerModels.parse_file(caseFile)
network_data = make_basic_network(network_data) # no dc lines, switches, inactive components, ...

# print_summary(network_data)
N_bus = length(network_data["bus"])
N_load = length(network_data["load"])
N_gen = length(network_data["gen"])
N_branch = length(network_data["branch"])

incidence = calc_basic_incidence_matrix(network_data) * -1 # we use inverse notation for directions
Zbr = calc_basic_branch_series_impedance(network_data)
# branch admittance
Ybr = Diagonal(inv.(Zbr))

# P and Q values for bus injection of all buses (these are not nec. PQ buses or contain loads)
# injection = calc_basic_bus_injection(network_data)
# P = real(injection)
# Q = imag(injection)

# P and Q values for load buses
# load_idx = getLoadIndices(network_data)
Pload, Qload = getLoadPQ(network_data)

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


# Initialization for 1 uncertainty (bus 8/load 5)
function initUncertainty_1(sys::Dict)
    α = 2
    β = 4.666
    op = Beta01OrthoPoly(maxDeg, α, β)

    println("Polynomial basis of size $(op.deg+1):")
    showbasis(op, digits = 2)
    println()

    # PCE of demand. Compute affine coefficients for each univariate unvertainty dimension and combine them
    pd = zeros(sys[:Nd], maxDeg + 1)
    qd = zeros(sys[:Nd], maxDeg + 1)
    pd[:, 1] = sys[:Pd]
    qd[:, 1] = sys[:Qd]
    pd[5, [1, 2]] = calculateAffinePCE(op) # random load 
    qd[5, [1, 2]] = copy(pd[5, [1,2]])

    global unc = Dict(:opq => op,
        :dim => op.deg + 1,
        :pd => pd,
        :qd => qd)

    global samples = sampleMeasure(numSamples, unc[:opq])  # NI model samples
    global ξ = sampleMeasure(5000, unc[:opq])  # Evaluation samples
end

# Initialization for 3 uncertainties
function initUncertainty_3()
    α = [3, 2]
    β = [2, 4]
    op1 = GaussOrthoPoly(maxDeg)
    op2 = Beta01OrthoPoly(maxDeg, α[1], β[1])
    op3 = Beta01OrthoPoly(maxDeg, α[2], β[2])

    global unc = setupUncertaintyMulti(maxDeg, [op1, op2, op3], 3)

    global samples = sampleMeasure(numSamples, unc[:opq])  # NI model samples
    global ξ = sampleMeasure(5000, unc[:opq])  # Evaluation samples
end


# Initialization for 4 uncertainties
function initUncertainty_4(numUnc::Int)
    # ...
end
