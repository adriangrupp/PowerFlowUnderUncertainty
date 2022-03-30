export showData,
compareCoefficients,
compareToMCMoments

### Evaluation of experiment results ###

"""
Show method to simply access the stored data
"""
function showData(file::String)
    f = load(file)
    haskey(f, "pf_state") && display(f["pf_state"])
    haskey(f, "moments")  && display(f["moments"])
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
    
    display(momentsMC)
    display(momentsPCE)

    for (key, val) in momentsMC
        if haskey(momentsPCE, key)
            mat1, mat2 = val, momentsPCE[key] # coefficients are stored as matrices
            diff = round.(mat1 - mat2, digits = 8) # difference between all entries
            # compare diff rowwise, i.e. for each bus
            for (i, row) in enumerate(eachrow(diff))
                println("($key, $i)  \t Error mean:\t", row[1], "\t\tError std:\t", row[2])
            end
        end
    end
end
