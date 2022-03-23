# Compare PCE coefficients
# compareCoefficients("coefficients/SPF_intrusive.jld", "coefficients/SPF_sparse.jld")
# compareCoefficients("coefficients/SPF_intrusive.jld", "coefficients/SPF_NI.jld")
# compareCoefficients("coefficients/SPF_NI_2unc.jld", "coefficients/SPF_sparse_2unc.jld")

# Compare Moments with MC
compareToMCMoments("coefficients/SPF_MC_moments.jld", "coefficients/SPF_NI_moments.jld")
compareToMCMoments("coefficients/SPF_MC_moments.jld", "coefficients/SPF_sparse_moments.jld")