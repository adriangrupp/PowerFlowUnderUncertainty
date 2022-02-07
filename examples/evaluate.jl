using JLD

### Evaluation of experiment results ###

# Comparison of PCE coefficients of different methods by ∞-Norm
function compareCoefficients(file1::String, file2::String)
    f1 = load(file1)
    f2 = load(file2)

    println("Comparing PCE coefficients between $(basename(file1)) and $(basename(file2))")

    coeffs1 = f1["pf_state"]
    coeffs2 = f2["pf_state"]

    # Loop over all PCE coefficient sets 
    for (key, val) in coeffs1
        if haskey(coeffs2, key)
            mat1, mat2 = val, coeffs2[key] # coefficients are stored as matrices

            diff = mat1 - mat2 # difference between coefficients
            # compare coefficients rowwise, i.e. for each bus
            for (i, row) in enumerate(eachrow(diff))
                println("∞-Norm for $key, $i:\t ", norm(diff, Inf))
            
            end
        end
    end
end

compareCoefficients("coefficients/SPF_intrusive.jld", "coefficients/SPF_sparse.jld")
# compareCoefficients("coefficients/SPF_intrusive.jld", "coefficients/SPF_NI.jld")