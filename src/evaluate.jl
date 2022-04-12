export showData,
    compareCoefficients,
    compareToMCMoments,
    numSampVsError

### Evaluation of experiment results ###

"""
Show method to simply access the stored data
"""
function showData(file::String)
    f = load(file)
    if haskey(f, "pf_state")
        JLD.display(f["pf_state"])
    elseif haskey(f, "moments")
        JLD.display(f["moments"])
    else
        display(f)
    end
end

"""
Comparison of PCE coefficients of different methods by ∞-Norm.
"""
function compareCoefficients(file1::String, file2::String)
    f1 = load(file1)
    f2 = load(file2)

    println("\nComparing PCE coefficients between $(basename(file1)) and $(basename(file2)):\n")

    coeffs1 = f1["pf_state"]
    coeffs2 = f2["pf_state"]

    # Loop over all PCE coefficient sets 
    for (key, val) in coeffs1
        if haskey(coeffs2, key)
            mat1, mat2 = val, coeffs2[key] # coefficients are stored as matrices

            # display(mat1)
            # display(mat2)

            diff = mat1 - mat2 # difference between coefficients
            # compare coefficients rowwise, i.e. for each bus
            for (i, row) in enumerate(eachrow(diff))
                println("∞-Norm for $(key), $(i):\t ", norm(row, Inf))

            end
        end
    end
end

"""
Compare PCE moments to Monte Carlo moments.
"""
function compareToMCMoments(mcFile::String, pceFile::String)
    f1 = load(mcFile)
    f2 = load(pceFile)

    println("\nComparing moments between $(basename(mcFile)) and $(basename(pceFile)):\n")

    momentsMC = f1["moments"]
    momentsPCE = f2["moments"]

    # display(momentsMC)
    # display(momentsPCE)

    for (key, val) in momentsMC
        if haskey(momentsPCE, key) && key == :e
            mat1, mat2 = val, momentsPCE[key] # Moments are stored as matrices
            diff = round.(mat1 - mat2, digits=10) # difference between all entries
            # compare diff rowwise, i.e. for each bus
            for (i, row) in enumerate(eachrow(diff))
                println("($key, $i)  \t Error mean:\t", row[1], "\t\tError std:\t", row[2])
            end
        end
    end
end

"""
Compare errors of all parameters vs numSamples in a plot. (moments or MSE)
"""
function numSampVsError(mcFile::String, pceFile1::String, pceFile2::String)
    f1 = load(mcFile)
    f2 = load(pceFile1)
    f3 = load(pceFile2)

    momentsMC = f1["moments"]

    ## Parse PCE Data
    # sort for numSamples (keys of dict)
    dictPCENi = sort(collect(f2), by=x -> parse(Int, x[1]))
    numSamp = [parse(Int, x[1]) for x in dictPCENi]
    dictPCESparse = sort(collect(f3), by=x -> parse(Int, x[1]))
    # Collection of all errors for all parameters(dict) -> all busses x numSamples(Matrix) -> Matrix entries are vectors of [numSamples, errMean, errStd]
    errorsNi = Dict{Symbol,Matrix{Vector}}()
    errorsSparse = Dict{Symbol,Matrix{Vector}}()
    m = length(numSamp)
    for (key, val) in momentsMC
        n = size(val, 1)
        errorsNi[key] = Matrix(undef, m, n)
        errorsSparse[key] = Matrix(undef, m, n)
    end

    ## Iterate PCE moments for each sample size (non-intrusive)
    for (i, momentsPCE) in enumerate(dictPCENi)
        momentsPCE = momentsPCE[2] # sub dictionaries store parameter data
        # Iterate momentsMC for each parameter
        for (key, val) in momentsMC
            # If parameter is also in PCE data, compute error difference
            if haskey(momentsPCE, key)
                mat1, mat2 = val, momentsPCE[key] # Moments are stored as matrices
                diff = mat1 - mat2
                diffVec = [diff[i, :] for i in 1:size(diff, 1)] # convert diff matrix to a vector of vectors
                # Store diff rowwise, i.e. for each bus together with the number of used samples
                errorsNi[key][i, :] = diffVec
            end
        end
    end

    ## Iterate PCE moments for each sample size (sparse)
    for (i, momentsPCE) in enumerate(dictPCESparse)
        momentsPCE = momentsPCE[2] # sub dictionaries store parameter data
        # Iterate momentsMC for each parameter
        for (key, val) in momentsMC
            # If parameter is also in PCE data, compute error difference
            if haskey(momentsPCE, key)
                mat1, mat2 = val, momentsPCE[key] # Moments are stored as matrices
                diff = mat1 - mat2
                diffVec = [diff[i, :] for i in 1:size(diff, 1)] # convert diff matrix to a vector of vectors
                # Store diff rowwise, i.e. for each bus together with the number of used samples
                errorsSparse[key][i, :] = diffVec
            end
        end
    end

    ## Plot diefference for all existing parameters for NI and sparse at once
    plotDir = "plots/10u_samples-moments_errors"
    for (key, errorMatNi) in errorsNi
        if haskey(errorsSparse, key)
            errorMatSparse = errorsSparse[key]
            if isassigned(errorMatNi, 1, 1) && isassigned(errorMatSparse, 1, 1) # if unassigned, don't plot
                plotSampleVsError(numSamp, errorMatNi, errorMatSparse, String(key), "moments", plotDir)
            end
        end
    end
end