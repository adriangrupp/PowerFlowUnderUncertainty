export getPQResult,
    getVoltageResult

"""
Return pg and qg for generator buses
"""
function getPQResult(pf)
    gen = sort(collect(pf["solution"]["gen"]), by = x -> parse(Int, x[1])) # sort solution/gen entries by IDs which are converted from string to Int
    pg = [el[2]["pg"] for el in gen]  # dicts are transferred into pairs. access values accodingly
    qg = [el[2]["qg"] for el in gen]
    return pg, qg
end
"""
Return v_re and v_im for generator buses
"""
function getVoltageResult(pf)
    bus = sort(collect(pf["solution"]["bus"]), by = x -> parse(Int, x[1]))
    vr = [el[2]["vr"] for el in bus]
    vi = [el[2]["vi"] for el in bus]
    return vr, vi
end