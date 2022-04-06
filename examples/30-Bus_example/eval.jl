using PowerFlowUnderUncertainty, JLD

# Compare PCE coefficients
# 1u: 
compareCoefficients("coefficients/SPF_NI.jld", "coefficients/SPF_sparse.jld")
# 2u:
compareCoefficients("coefficients/SPF_2u_NI.jld", "coefficients/SPF_2u_sparse.jld")
# 10u:
compareCoefficients("coefficients/SPF_10u_NI.jld", "coefficients/SPF_10u_sparse.jld")


# Compare Moments with MC
# 1u:
compareToMCMoments("coefficients/SPF_MC_moments.jld", "coefficients/SPF_NI_moments.jld")
compareToMCMoments("coefficients/SPF_MC_moments.jld", "coefficients/SPF_sparse_moments.jld")
# 2u:
compareToMCMoments("coefficients/SPF_2u_MC_moments.jld", "coefficients/SPF_2u_NI_moments.jld")
compareToMCMoments("coefficients/SPF_2u_MC_moments.jld", "coefficients/SPF_2u_sparse_moments.jld")
# 10u:
compareToMCMoments("coefficients/SPF_10u_MC_moments.jld", "coefficients/SPF_10u_NI_moments.jld")
compareToMCMoments("coefficients/SPF_10u_MC_moments.jld", "coefficients/SPF_10u_sparse_moments.jld")