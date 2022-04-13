using PowerFlowUnderUncertainty

# Compare PCE coefficients
# 1u: 
compareCoefficients("coefficients/SPF_NI.jld2", "coefficients/SPF_sparse.jld2")
# 2u:
compareCoefficients("coefficients/SPF_2u_NI.jld2", "coefficients/SPF_2u_sparse.jld2")
# 10u:
compareCoefficients("coefficients/SPF_10u_NI.jld2", "coefficients/SPF_10u_sparse.jld2")


# Compare Moments with MC
# 1u:
compareToMCMoments("coefficients/SPF_MC_moments.jld2", "coefficients/SPF_NI_moments.jld2")
compareToMCMoments("coefficients/SPF_MC_moments.jld2", "coefficients/SPF_sparse_moments.jld2")
# 2u:
compareToMCMoments("coefficients/SPF_2u_MC_moments.jld2", "coefficients/SPF_2u_NI_moments.jld2")
compareToMCMoments("coefficients/SPF_2u_MC_moments.jld2", "coefficients/SPF_2u_sparse_moments.jld2")
# 10u:
compareToMCMoments("coefficients/SPF_10u_MC_moments.jld2", "coefficients/SPF_10u_NI_moments.jld2")
compareToMCMoments("coefficients/SPF_10u_MC_moments.jld2", "coefficients/SPF_10u_sparse_moments.jld2")


# Compare iterated PCE methods
# 10u error wrt MC moments
numSampVsError("coefficients/SPF_10u_MC_moments.jld2", "coefficients/SPF_10u_NI_moments-iter.jld2", "coefficients/SPF_10u_sparse_moments-iter.jld2")
# 10u error wrt MSE
numSampVsError("coefficients/SPF_10u_NI_mse-iter.jld2", "coefficients/SPF_10u_sparse_mse-iter.jld2")