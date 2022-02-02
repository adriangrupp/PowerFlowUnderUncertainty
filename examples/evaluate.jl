using JLD

### Evaluation of experiment results ###

function compareCoefficients(file1::String, file2::String)
    f1 = load(file1)
    f2 = load(file2)

    println(f1)
end


compareCoefficients("coefficients/SPF_NI.jld", "coefficients/SPF_NI.jld")