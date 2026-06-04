# ==============================================================================
# CR_GOF_Simulation_Scenario6.r
#
# Simulation scenario 6: Models where heterogeneity is modeled with a covariate 
# and a random effect.
#
# We evaluate model GOFs (for 2 models only) using DHARMa. GOF residual plots 
# are assessed visually and we use the KS, Dispersion and Outlier tests. 
#
# We also report some summaries for population size estimates for each model 
# using partial likelihood with a Horvitz-Thompson estimator. We also include 
# population size estimates fitting a (Huggins) conditional 
# likelihood/Horvitz-Thompson estimator M_bh (GLM) using TMB.
#
# A full description of these data and their analysis appears in the 
# Supplementary Materials (Web Section C6) of the main text.
# ==============================================================================

library(DHARMa)
library(TMB)
library(tidyverse)
library(glmmTMB)

source("CR_GOF_Functions.R")

# The code below is required for fitting the M_h GLMM. It is a fully marginalized 
# CL approach that uses Laplace approximation via TMB for parameter estimation.

filename_XRE <- "TMB_Mh_RE.cpp"  # Compile C++ file.
modelname_XRE <- "TMB_Mh_RE"
TMB::compile(filename_XRE, 
             flags = "-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign")

dyn.load(dynlib(modelname_XRE))

nsim <- 50   # Number of simulations.

# Set parameters.

N <- 250      # True population size.
tau <- 10       # Number of capture occasions.

#beta_0 <- -0.5
#beta_1 <- 0.75
#sigma <- 1.25   # Standard deviation for the random effect.

beta_0 <- -1
beta_1 <- 1.3

sigma <- 1   # Standard deviation for the random effect.

res <- matrix(NA, nsim, 4)

folder_name <- "Sim6_GOF"
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}

set.seed(123)

for (ii in 1:nsim) {

  X <- runif(N, -2, 2)  # Covariate to always be included in the model.
  
  Z <- rnorm(N, 0, sd = sigma)  # Random effects (epsilon_i) included in the model.
  
  # Capture probability.
  
  P_i <- plogis(beta_0 + beta_1*X + Z) # Linear with random effect.
  
  # Generate a single capture-recapture data set here.
  
  cap.hist0 <- matrix(rbinom(N*tau, 1, P_i), nrow = N, ncol = tau)
  
  cap.hist <- cap.hist0[rowSums(cap.hist0) != 0, ] # Remove zero (uncaptured) rows.
  
  D <- nrow(cap.hist)               # Unique number of individuals captured.
  
  X_i <- X[rowSums(cap.hist0) != 0]     # Covariate for captured individuals.
  
  x.obs <- as.matrix(X_i)
  colnames(x.obs) <- c("x.obs")
  
  # Construct (Binomial) recapture-based quantities under the partial likelihood 
  # for fitting M_h-only models.
  
  pl_data <- prepare_recap_data_Mh(cap.hist, x.obs, tau)
  
  # Start fitting M_h models here.
  
  # Incorrectly specified: GLM with covariate but without the random effect.
  
  model1 <- glm(cbind(Y_recap, n_opp - Y_recap) ~ x.obs,
                family = binomial(link = "logit"), data = pl_data)
  
  residuals_output1 <- simulateResiduals(model1, n = 1000)  # DHARMa simulated residuals.
  
  # Get population size estimates.
  
  p_i_hat <- predict(model1, newdata = data.frame(x.obs = X_i), type = "response")
  pi_i <- 1 - (1 - p_i_hat)^tau
  N_hat1 <- sum(1/pi_i)     # The Horvitz-Thompson estimator.
  
  # Also incorrectly specified: A GLMM with individual random effects (and covariate).
  
  model2 <- glmmTMB(cbind(Y_recap, n_opp - Y_recap) ~ x.obs + (1 | ID), 
                  family = binomial(link = "logit"), data = pl_data)
  
  residuals_output2 <- simulateResiduals(model2, n = 1000)  # DHARMa simulated residuals.
  
  # Get population size estimates.
  
  fix_ef <- fixef(model2)
  
  # Calculate linear predictor for ALL D individuals (using X_i). Start with the 
  # fixed effect portion: theta0 + theta1*X_i.
  
  eta_i <- fix_ef[1]$cond[1] + fix_ef[1]$cond[2]*X_i
  
  # Add the Random Effects (intercept shifts) for those we have them for. For 
  # individuals caught before the last trial.
  
  t1 <- apply(cap.hist, 1, function(x) match(1, x)) 
  caught_last <- (t1 == tau)
  
  re_modes <- ranef(model2)
  eta_i[!caught_last] <- eta_i[!caught_last] + c(re_modes$cond$ID$`(Intercept)`)
  p_i_hat <- plogis(eta_i)
  pi_i <- 1 - (1 - p_i_hat)^tau
  
  N_hat2 <- sum(1/pi_i)  # The Horvitz-Thompson estimator.
  
  # A full marginalized approach using the conditional likelihood.
  
  X_design <- cbind(1, X_i) 
  
  tidbits_data <- list(
    y = as.matrix(cap.hist), 
    X = as.matrix(X_design)
  )
  
  tidbits_parameters <- list(
    betas = c(fix_ef[1]$cond[1], fix_ef[1]$cond[2]),
    individual_effects = rep(0, D),
    logsigma_individual_effects = 0 
  )
  
  tidbits_random <- c("individual_effects")
  
  objs <- MakeADFun(data = tidbits_data, 
                    parameters = tidbits_parameters, 
                    random = tidbits_random, 
                    DLL = modelname_XRE, silent = TRUE)
  
  fit_XRE_marg <- nlminb(start = objs$par, objective = objs$fn, gradient = objs$gr, 
                           control = list(trace = 0, iter.max = 5000, eval.max = 5000))
  
  fit_XRE_marg_results <- sdreport(objs, bias.correct = FALSE, 
                                     ignore.parm.uncertainty = FALSE, skip.delta.method = FALSE)
  fit_XRE_marg_estimates <- as.list(fit_XRE_marg_results, what = "Estimate", report = TRUE)
  fit_XRE_marg_all_results <- summary(fit_XRE_marg_results, 
                                        select = "report", p.value = TRUE) %>% as.data.frame() %>% rownames_to_column(var = "Parameters")
  
  rm(fit_XRE_marg, fit_XRE_marg_results)
  
  # Get population size.
  
  N_hat3 <- fit_XRE_marg_all_results %>% dplyr::filter(Parameters == "population_size")
  
  res[ii, ] <- c(D, N_hat1, N_hat2, N_hat3[2]$`Estimate`)
  
  # DHARMa plots: Side-by-side plots.
  
  figname <- file.path(folder_name, paste0("sim6_GOF_plots_dataset_", ii, ".png"))
  
  png(file = figname, width = 11, height = 9, units = "in", res = 300)
  
  par(mfrow = c(2, 2))
  
  plotQQunif(residuals_output1, main = paste("Dataset_", ii, ": Model 1 (GLM)"))
  plotResiduals(residuals_output1)
  plotQQunif(residuals_output2, main = paste("Dataset_", ii, ": Model 2 (GLMM)"))
  plotResiduals(residuals_output2)
  
  dev.off()
}

#...............................................................................

res_all <- cbind(apply(res, 2, mean), apply(res, 2, median), 
                 apply(res, 2, sd), apply(res, 2, mad), 
                 sqrt(apply((res - N)^2, 2, mean)))

colnames(res_all) <- c("Mean(N.hat)", "Median(N.hat)", 
                       "SD(N.hat)", "MAD(N.hat)", "RMSE")
rownames(res_all) <- c("D", "PL-GLM", "PL-GLMM", "CL-GLMM (marg)")
res_all

# Relative % bias for N:

((res_all[, 2] - N)/N)*100
