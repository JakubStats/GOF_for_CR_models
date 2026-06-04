# ==============================================================================
# CR_GOF_Simulation_Scenario4.r
#
# Simulation scenario 4: Models with time/capture effects and heterogeneity.
#
# We evaluate model GOFs (for 2 models only) using DHARMa. GOF residual plots 
# are assessed visually and we use the KS and Dispersion. 
#
# We also report some summaries for population size estimates for each model 
# using partial likelihood with a Horvitz-Thompson estimator. We also include 
# population size estimates fitting a (Huggins) conditional 
# likelihood/Horvitz-Thompson estimator M_th (GLM) using VGAM.
#
# A full description of these data and their analysis appears in the 
# Supplementary Materials (Web Section C4).
# ==============================================================================

library(DHARMa)
library(VGAM)

source("CR_GOF_Functions.R")

nsim <- 50   # Number of simulations.

# Set parameters.

N <- 50      # True population size.
tau <- 15       # Number of capture occasions.

if(tau != 7 && tau != 15) stop("Unsupported number of capture occasions (tau). 
                               Tau needs to be either 7 or 15 for these sims.")

if (tau == 7) alpha <- c(1, 3, 1, -3, 0, 2, -2)  # Time effects.
if (tau == 15) alpha <- c(1, 3, 1, -3, 0, 2, -2, -2, 0, -3, 1, 0, 3, -1, 2, 1, 0)  # Time effects.

beta_0 <- -2.2   # Intercept.
beta_1 <- -1.5   # Covariate coef.

res <- matrix(NA, nsim, 4)

folder_name <- "Sim4_GOF"
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}

set.seed(123)

for (ii in 1:nsim) {

  X <- runif(N, -2, 2)  # A individual covariate added as heterogeneity.
  
  pr <- matrix(0, ncol = tau, nrow = N)
  y1 <- matrix(0, ncol = tau + 1, nrow = N)
  
  for (i in 1:N) {
    for (j in 1:tau) {
      pr[i, j] <- plogis(beta_0 + alpha[j] + beta_1*X[i])
      
      y1[i, j] <- rbinom(1, size = 1, prob = pr[i, j])
    }
    
    y1[i, tau + 1] <- sum(y1[i, ])
  }
  
  cap.hist <- y1[y1[, (tau + 1)] != 0, 1:tau] # Filter only captured, remove sum column.
  X_i <- X[y1[, (tau + 1)] != 0]  # Filter covariates for those same individuals.
  D <- nrow(cap.hist)           
  
  x.obs <- as.matrix(X_i)
  colnames(x.obs) <- c("x.obs")
  
  # Construct (Bernoulli) recapture-based quantities under the partial likelihood
  # for fitting M_t- and M_b-type models.
  
  pl_data <- prepare_recap_data_Mtbh(as.matrix(cap.hist), x.obs, tau)
  
  # Incorrectly specified: GLM (M_h) using the covariate only.
  
  mod_mh <- glm(Y_recap ~ x.obs, family = binomial(link = "logit"), 
                data = pl_data)
  
  preds_by_id <- split(predict(mod_mh, type = "response"), pl_data$ID)
  pi_i <- sapply(preds_by_id, function(x) 1 - prod(1 - x))
  N_hat1 <- sum(1/pi_i)
  
  residuals_output1 <- simulateResiduals(mod_mh, n = 1000)
  res_grouped1 <- recalculateResiduals(residuals_output1, group = pl_data$ID)
  
  # Correctly specified: As above, but now account for time effects. That is,
  # fit a M_th GLM.
  
  mod_mth <- glm(Y_recap ~ x.obs + occasion, family = binomial(link = "logit"), 
                 data = pl_data)
  
  mth_grid <- expand.grid(
    ID = 1:D,
    occ_idx = factor(1:tau, levels = 2:tau)
  )
  
  mth_grid$x.obs <- X_i[mth_grid$ID]
  
  # Handle the 'occasion' factor levels carefully. Here 'occasion' starts from 
  # 2 (the first possible recapture). We must provide a value for Day 1. 
  # Since mod_mth uses occasion as a factor, Day 1 is considered the 
  # "Reference Level" (the intercept).
  
  mth_grid$occasion <- factor(mth_grid$occ_idx, levels = 2:tau)
  
  # Predict probabilities for every individual on every day
  # Note: For Day 1, predict() will use the intercept (the reference level).
  
  p_all_days <- predict(mod_mth, newdata = mth_grid, type = "response")
  p_all_days[is.na(p_all_days)] <- plogis(coef(mod_mth)[1] + 
                                            coef(mod_mth)["x.obs"]*mth_grid$x.obs[is.na(p_all_days)])
  
  # Group by ID and calculate the full Huggins Pi (at least once in tau days).
  
  p_list <- split(p_all_days, mth_grid$ID)
  pi_i <- sapply(p_list, function(x) 1 - prod(1 - x, na.rm = TRUE))
  
  N_hat2 <- sum(1/pi_i)
  
  residuals_output2 <- simulateResiduals(mod_mth, n = 1000)
  res_grouped2 <- recalculateResiduals(residuals_output2, group = pl_data$ID)
  
  # Fit model M_th with VGAM (which uses the Huggins' conditional likelihood) 
  # to get the pop size.
  
  data_Mth_VGLM <- data.frame(cbind(cap.hist, X_i))
  
  if (tau == 7) {
    colnames(data_Mth_VGLM) <- c("y1", "y2", "y3", "y4", "y5", "y6", "y7", "X_i")
    
    M_th <- vglm(cbind(y1, y2, y3, y4, y5, y6, y7) ~ X_i, 
                 posbernoulli.t, data = data_Mth_VGLM, trace = FALSE) 
  }
  
  if (tau == 15) {
    colnames(data_Mth_VGLM) <- c("y1", "y2", "y3", "y4", "y5", "y6", "y7", "y8",
                                 "y9", "y10", "y11", "y12", "y13", "y14", "y15", "X_i")
    
    M_th <- vglm(cbind(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, 
                       y14, y15) ~ X_i, 
                 posbernoulli.t, data = data_Mth_VGLM, trace = FALSE) 
  }
  
  N_hat3 <- M_th@extra$N.hat
  
  # DHARMa plots: Side-by-side plots.
  
  figname <- file.path(folder_name, paste0("sim4_GOF_plots_dataset_", ii, ".png"))
  
  png(file = figname, width = 11, height = 9, units = "in", res = 300)
  
  par(mfrow = c(2, 2))
  
  plotQQunif(res_grouped1, 
             main = paste("Dataset_", ii, ": Model 1 (Misspecified, M_h)"))
  plotResiduals(res_grouped1)
  plotQQunif(res_grouped2, 
             main = paste("Iter", ii, ": Model 2 (Correctly specified, M_th)"))
  plotResiduals(res_grouped2)
  
  dev.off()
  
  res[ii, ] <- c(D, N_hat1, N_hat2, N_hat3)
}

#...............................................................................

res_all <- cbind(apply(res, 2, mean), apply(res, 2, median), 
                 apply(res, 2, sd), apply(res, 2, mad), 
                 sqrt(apply((res - N)^2, 2, mean)))

colnames(res_all) <- c("Mean(N.hat)", "Median(N.hat)", 
                       "SD(N.hat)", "MAD(N.hat)", "RMSE")
rownames(res_all) <- c("D", "PL-M_h", "PL-M_th", "CL-M_th")
res_all

# Relative % bias for N:

((res_all[, 2] - N)/N)*100
