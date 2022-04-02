export runPfModel,
runOpfModel


"""
Read in network data
"""
function readCaseFlie(caseFile)
    PowerModels.logger_config!("error") # Suppress warnings on phase angles
    network_data = PowerModels.parse_file(caseFile)
    network_data = make_basic_network(network_data) # basuc network: no dc lines, switches, inactive components, ...
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
Managing of the PF model. Wrapper function for NI-algo.
Input:  network_data - modified netword description of PQ values. solver - initialized NLP solver
Return: dict of PF outputs.
"""
function runPfModel(network_data, solver)
    pf = PowerModels.run_pf(network_data, ACRPowerModel, solver)

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

"""
Managing of the OPF model. Wrapper function for NI-algo.
Input:  x - sampled value for active power of PQ bus.
Return: dict of OPF outputs.
"""
function runOpfModel(network_data, solver)
    pf = PowerModels.run_opf(network_data, ACRPowerModel, solver)

    # check if solution was feasible
    status = pf["termination_status"]
    status == OPTIMAL || status == LOCALLY_SOLVED ? nothing : error("Potentially no solution found: ", status)

    pg, qg = getGenPQResult(pf)
    vr, vi = getVoltageResult(pf)

    objective = pf["objective"] # objective function value

    return Dict(:pg => pg,
        :qg => qg,
        :e => vr, # we use e and f as convention for real and imaginary voltage parts
        :f => vi,
        :obj => objective)
end