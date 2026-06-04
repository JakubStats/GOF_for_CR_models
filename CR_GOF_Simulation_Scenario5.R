# ==============================================================================
# CR_GOF_Simulation_Scenario5.r
#
# Simulation scenario 5: Models with behavioural effects, time/capture effects 
# and heterogeneity.
#
# We evaluate model GOFs (for 3 models) using DHARMa. GOF residual plots 
# are assessed visually and we use the KS and Dispersion. 
#
# We also report some summaries for population size estimates for each model 
# using partial likelihood with a Horvitz-Thompson estimator. We also include 
# population size estimates fitting a (Huggins) conditional 
# likelihood/Horvitz-Thompson estimator M_tbh (GLM) using VGAM.
#
# A full description of these data and their analysis appears in the 
# Supplementary Materials (Web Section C5).
# ==============================================================================

library(DHARMa)
library(VGAM)

source("CR_GOF_Functions.R")

nsim <- 50   # Number of simulations.

# Set parameters.

N <- 50        # True population size.
tau <- 7      # Number of capture occasions.

if(tau != 7 && tau != 15) stop("Unsupported number of capture occasions (tau). 
                               Tau needs to be either 7 or 15 for these sims.")

if (tau == 7) alpha <- c(1, 3, 1, -3, 0, 2, -2)  # Time effects.
if (tau == 15) alpha <- c(1, 3, 1, -3, 0, 2, -2, -2, 0, -3, 1, 0, 3, -1, 2, 1, 0)  # Time effects.

beta_0 <- -2.2   # Intercept.
beta_1 <- -1.5   # Covariate coef.

rho_val <- -2          # Behave. effect (trap happy +ve and trap shy -ve).
rho <- rep(rho_val, N)

res <- matrix(NA, nsim, 5)

folder_name <- "Sim5_GOF"
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}

set.seed(123)

for (ii in 1:nsim) {

  X <- runif(N, -2, 2)  # A individual covariate added as heterogeneity.
  
  eta0 <- beta_0 + X*beta_1
  eta1 <- outer(eta0, alpha, "+")
  eta2 <- outer((eta0 + rho), alpha, "+")
  pij1 <- plogis(eta1)
  pij2 <- plogis(eta2)             
  
  #ru <- matrix(runif(N*tau), ncol = tau)
  #hist1 <- matrix(as.numeric(ru < pij1), ncol = tau)
  #hist2 <- matrix(as.numeric(ru < pij2), ncol = tau)
  #histid <- matrix(as.numeric(t(apply(hist1, 1, cumsum)) > 0), ncol = tau)
  #tot.cap <- histid[, tau] > 0
  #histid <- cbind(rep(0, N), histid[, 1:(tau - 1)])
  #cc <- hist1*(1 - histid) + hist2*histid
  
  # First occasion is always based on the 'initial' probability (pij1). This
  # sets up the Markov behavioral response to capture.
  
  cc <- matrix(0, nrow = N, ncol = tau)
  cc[, 1] <- as.numeric(runif(N) < pij1[, 1])
  
  # Subsequent occasions depend on the outcome of the previous one.
  
  # If caught at j-1, use pij2 (the trap response probability).
  # If NOT caught at j-1, use pij1 (the standard probability).
  
  for (j in 2:tau) {
    prob_current <- ifelse(cc[, j - 1] == 1, pij2[, j], pij1[, j])
    cc[, j] <- as.numeric(runif(N) < prob_current)
  }
  
  tot.cap <- rowSums(cc) > 0
  cap.hist <- cc[tot.cap, ]      # Observed capture history matrix.
  
  # Generate a single capture-recapture data set here.
  
  D <- nrow(cap.hist)   # Number of (unique) individuals captured.
  
  X_i <- X[tot.cap]   # Observed covariate for captured individuals.    
  
  x.obs <- as.matrix(X_i)
  colnames(x.obs) <- c("x.obs")
  
  # Construct (Bernoulli) recapture-based quantities under the partial likelihood
  # for fitting M_t- and M_b-type models.
  
  pl_data <- prepare_recap_data_Mtbh(as.matrix(cap.hist), x.obs, tau)
  
  # Incorrectly specified: Fit a M_bh-type GLM (accounts for behav. response and hete.).
  
  mod_mbh <- glm(Y_recap ~ x.obs + cap_prev, family = binomial(link = "logit"), 
                 data = pl_data)
  
  mbh_study_grid <- expand.grid(
    ID = 1:D,
    occ = 1:tau  # Ensure it projects over all occasions.
  )
  
  # For Huggins M_bh, we assume cap_prev = 0 for the "potential" probability of 
  # the first capture.
  
  mbh_study_grid$x.obs<- X_i[mbh_study_grid$ID] # Ensure X_i matches the ID.
  mbh_study_grid$cap_prev <- 0 # Standard Huggins assumes 0 for the PI calculation.
  
  # Predict across the whole study.
  
  mbh_preds <- predict(mod_mbh, newdata = mbh_study_grid, type = "response")
  mbh_preds_list <- split(mbh_preds, mbh_study_grid$ID)
  
  pi_i <- sapply(mbh_preds_list, function(x) 1 - prod(1 - x))
  N_hat1 <- sum(1/pi_i)           # The Horvitz-Thompson estimator.
  
  residuals_output1 <- simulateResiduals(mod_mbh, n = 1000)
  res_grouped1 <- recalculateResiduals(residuals_output1, group = pl_data$ID)
  
  # Incorrectly specified: Fit a M_th-type GLM (accounts for time effects and hete.).
  
  mod_mth <- glm(Y_recap ~ x.obs + occasion, family = binomial(link = "logit"), 
                 data = pl_data)
  
  mth_grid <- expand.grid(
    ID = 1:D,
    occ_idx = factor(1:tau, levels = 2:tau)
  )
  
  mth_grid$x.obs <- X_i[mth_grid$ID]
  mth_grid$occasion <- factor(mth_grid$occ_idx, levels = 2:tau)
  
  p_all_days <- predict(mod_mth, newdata = mth_grid, type = "response")
  p_all_days[is.na(p_all_days)] <- plogis(coef(mod_mth)[1] + 
                                            coef(mod_mth)["x.obs"]*mth_grid$x.obs[is.na(p_all_days)])
  p_list <- split(p_all_days, mth_grid$ID)
  pi_i <- sapply(p_list, function(x) 1 - prod(1 - x, na.rm = TRUE))
  
  N_hat2 <- sum(1/pi_i)      # The Horvitz-Thompson estimator.
  
  residuals_output2 <- simulateResiduals(mod_mth, n = 1000)
  res_grouped2 <- recalculateResiduals(residuals_output2, group = pl_data$ID)
  
  # Correctly specified: Fit the full M_tbh model (account for all three sources).
  
  mod_mtbh <- glm(Y_recap ~ x.obs + cap_prev + occasion, 
                  data = pl_data, family = binomial)
  
  mtbh_grid <- expand.grid(
    ID = 1:D,
    occ_idx = factor(1:tau, levels = 2:tau)
  )
  
  mtbh_grid$x.obs <- X_i[mtbh_grid$ID]
  mtbh_grid$occasion <- factor(mtbh_grid$occ_idx, levels = 2:tau)
  mtbh_grid$cap_prev <- 0
  
  p_all_days_mtbh <- predict(mod_mtbh, newdata = mtbh_grid, type = "response")
  p_all_days_mtbh[is.na(p_all_days_mtbh)] <- plogis(coef(mod_mtbh)[1] + 
                                                      coef(mod_mtbh)["x.obs"]*mtbh_grid$x.obs[is.na(p_all_days_mtbh)])
  p_list_mtbh <- split(p_all_days_mtbh, mtbh_grid$ID)
  pi_i_mtbh <- sapply(p_list_mtbh, function(x) 1 - prod(1 - x, na.rm = TRUE))
  
  N_hat3 <- sum(1/pi_i_mtbh)     # The Horvitz-Thompson estimator.
  
  residuals_output3 <- simulateResiduals(mod_mtbh, n = 1000)
  res_grouped3 <- recalculateResiduals(residuals_output3, group = pl_data$ID)
  
  # Fit model M_tbh with VGAM (which uses the Huggins' conditional likelihood) 
  # to get the pop size.
  
  data_Mtbh_VGLM <- data.frame(cbind(cap.hist, X_i))
  
  if (tau == 7) {
    colnames(data_Mtbh_VGLM) <- c("y1", "y2", "y3", "y4", "y5", "y6", "y7", "X_i")
  
  M_tbh <- vglm(cbind(y1, y2, y3, y4, y5, y6, y7) ~ X_i, 
               posbernoulli.tb, data = data_Mtbh_VGLM, trace = FALSE) 
  }
  
  if (tau == 15) {
    colnames(data_Mtbh_VGLM) <- c("y1", "y2", "y3", "y4", "y5", "y6", "y7", "y8",
                                  "y9", "y10", "y11", "y12", "y13", "y14", "y15", "X_i")
    
    M_tbh <- vglm(cbind(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, 
                        y13, y14, y15) ~ X_i, 
                  posbernoulli.tb, data = data_Mtbh_VGLM, trace = FALSE) 
    
  }
  
  N_hat4 <- M_tbh@extra$N.hat
  
  # DHARMa plots: Side-by-side plots.
  
  figname <- file.path(folder_name, paste0("sim5_GOF_plots_dataset_", ii, ".png"))
  
  png(file = figname, width = 11, height = 9, units = "in", res = 300)
  
  par(mfrow = c(3, 1))
  
  plotQQunif(res_grouped1, 
             main = paste("Dataset_", ii, ": Model 1 (Misspecified, M_bh)"))
  plotQQunif(res_grouped2, 
             main = paste("Dataset_", ii, ": Model 2 (Misspecified, M_th)"))
  plotQQunif(res_grouped3, 
             main = paste("Dataset_", ii, ": Model 3 (Correctly specified, M_tbh)"))
  
  dev.off()
  
  res[ii, ] <- c(D, N_hat1, N_hat2, N_hat3, N_hat4)
}

#...............................................................................

res_all <- cbind(apply(res, 2, mean), apply(res, 2, median), 
                 apply(res, 2, sd), apply(res, 2, mad), 
                 sqrt(apply((res - N)^2, 2, mean)))

colnames(res_all) <- c("Mean(N.hat)", "Median(N.hat)", "SD(N.hat)", "MAD(N.hat)", "RMSE")
rownames(res_all) <- c("D", "PL-M_bh", "PL-M_th", "PL-M_tbh", "CL-M_tbh")
res_all

# Relative % bias for N:

((res_all[, 2] - N)/N)*100
