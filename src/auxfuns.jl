export getLoadIndices,
getCost,
getLoadPQ,
getGenPQResult,
getVoltageResult


"""
Get indices of load buses
"""
function getLoadIndices(network_data)
    idx = [el[2]["load_bus"] for el in network_data["load"]]
    return sort(idx)
end
"""
Get quadradic and linear generator cost
"""
function getCost(network_data)
    costquad = [el[2]["cost"][1] for el in network_data["gen"]]
    costlin  = [el[2]["cost"][2] for el in network_data["gen"]]
    return costquad, costlin
end
"""
Get P and Q values for load buses
"""
function getLoadPQ(network_data)
    P = [el[2]["pd"] for el in network_data["load"]]
    Q = [el[2]["qd"] for el in network_data["load"]]
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