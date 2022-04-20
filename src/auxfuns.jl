export getLoadIndices,
getCost,
getLoadPQ,
getGenPQ,
getGenPQResult,
getVoltageResult


"""
Get indices of all load buses
"""
function getLoadIndices(network_data)
    load = sort(collect(network_data["load"]), by=x -> parse(Int, x[1])) # sort solution/gen entries by IDs which are converted from string to Int
    idx = [el[2]["load_bus"] for el in load]
    return sort(idx)
end
"""
Get quadradic and linear generator cost
"""
function getCost(network_data)
    gen = sort(collect(network_data["gen"]), by=x -> parse(Int, x[1])) # sort solution/gen entries by IDs which are converted from string to Int
    costquad = [el[2]["cost"][1] for el in gen]
    costlin = [el[2]["cost"][2] for el in gen]
    return costquad, costlin
end
"""
Get P and Q values for load buses
"""
function getLoadPQ(network_data)
    load = sort(collect(network_data["load"]), by = x -> parse(Int, x[1])) # sort solution/gen entries by IDs which are converted from string to Int
    P = [el[2]["pd"] for el in load]
    Q = [el[2]["qd"] for el in load]
    return P, Q
end
"""
Get P and Q values for generator buses
"""
function getGenPQ(network_data)
    gen = sort(collect(network_data["gen"]), by = x -> parse(Int, x[1])) # sort solution/gen entries by IDs which are converted from string to Int
    P = [el[2]["pg"] for el in gen]
    Q = [el[2]["qg"] for el in gen]
    return P, Q
end
"""
Get pg and qg for generator buses
"""
function getGenPQResult(pf)
    gen = sort(collect(pf["solution"]["gen"]), by = x -> parse(Int, x[1])) # sort solution/gen entries by IDs which are converted from string to Int
    pg = [el[2]["pg"] for el in gen]  # dicts are transferred into pairs. access values accodingly
    qg = [el[2]["qg"] for el in gen]
    return pg, qg
end
"""
Get v_re and v_im for all buses
"""
function getVoltageResult(pf)
    bus = sort(collect(pf["solution"]["bus"]), by = x -> parse(Int, x[1]))
    vr = [el[2]["vr"] for el in bus]
    vi = [el[2]["vi"] for el in bus]
    return vr, vi
end