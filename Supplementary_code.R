library(metafor)

#############################################
# Setup
#############################################

# In meta-analysis, researchers normally have available a set of observed effect sizes and corresponding variances

# We show that here for k = 5 studies
k <- 5

# Observed effect sizes
Ti <- c(0.1, 0.2, 0.3, 0.4, 0.5)

# Observed effect size variances
vi <- c(0.005, 0.002, 0.01, 0.009, 0.006)

# Sample size in each study
ni <- c(200, 500, 100, 111, 167)

# Population size of sample in each study
Ni <- c(4000, 10000, 1000, 3000, 2500)

#############################################
# Simple Random Sampling (SRS)
#############################################
# SRS weights
w_srs <- 1/vi

# SRS summary effect size estimate
T_srs <- sum((w_srs*Ti))/sum(w_srs)

# SRS summary effect size variance estimate
Var_T_srs <- 1 / sum(w_srs)

# Alternatively, as SRS uses inverse-variance weights, we can use standard CE (called FE in metafor) MA 
rma(yi=Ti, vi=vi, method="FE")

#############################################
# Stratified Simple Random Sampling (SSRS)
#############################################
# SSRS weights
w_ssrs <- Ni

# SSRS summary effect size estimate
T_ssrs <- sum((w_ssrs*Ti))/sum(w_ssrs)

# SSRS summary effect size variance estimate
Var_T_ssrs <- sum((w_ssrs^2)*vi) / (sum(w_ssrs)^2)

#############################################
# Cluster Sampling (CS)
#############################################

# CS weights
w_cs <- Ni

# SSRS summary effect size estimate
T_cs <- sum((w_cs*Ti))/sum(w_cs)

# Variance calculation requires the "true" number of clusters in the population
# This may be unknown or clusters may be hypothetical, so instead use RVE
# Here with k=5 studies, we use RVE with HC2 adjustment
Ai <- (1 - (w_cs / sum(w_cs)))^-1
Var_T_cs <- sum(Ai * (w_cs^2) * ((Ti - T_cs)^2)) / (sum(w_cs)^2)






