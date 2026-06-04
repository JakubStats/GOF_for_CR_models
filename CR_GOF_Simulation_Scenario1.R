# ==============================================================================
# CR_GOF_Simulation_Scenario1.r
#
# Simulation scenario 1: Models where heterogeneity is modeled with covariates.
#
# We evaluate model GOFs (for 2 models only) using DHARMa. GOF residual plots 
# are assessed visually and we use the KS, Dispersion and Outlier tests. 
#
# We repeat this 50 times. 
#
# We also report some summaries for population size estimates for each model 
# using partial likelihood with a Horvitz-Thompson estimator. We also include 
# population size estimates fitting a (Huggins) conditional 
# likelihood/Horvitz-Thompson estimator M_h (GLM) using VGAM.
#
# A full description of these data and their analysis appears in the 
# Supplementary Materials (Web Section C1).
# ==============================================================================

library(DHARMa)
library(VGAM)

source("CR_GOF_Functions.R")

nsim <- 50   # Number of simulations.

# Set parameters.

# For W ~ normal dist. use:

N <- 50       # True population size.
tau <- 7       # Number of capture occasions.

if(tau != 7 && tau != 15) stop("Unsupported number of capture occasions (tau). 
                               Tau needs to be either 7 or 15 for these sims.")

# IMPORTANT: Toggle below to change the distribution of W.

dist.W <- "Normal"
#dist.W <- "Binary"

# For W ~ normal dist. use:

if (dist.W == "Normal") {
  beta_0 <- 1
  beta_1 <- 1
  beta_2 <- -1
}

# For W ~ binomial dist. use:

if (dist.W == "Binary") {
  beta_0 <- -0.5
  beta_1 <- 1.5
  beta_2 <- -2
}

res <- matrix(NA, nsim, 4)

folder_name <- "Sim1_GOF"
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}

set.seed(123)

for (ii in 1:nsim) {
  
  # Covariate to always be included in the model.
  
  X <- rnorm(N, 0, 1)
  
  # Include an additional covariate.
  
  # IMPORTANT: Toggle here to change dist. of W.
  
  if (dist.W == "Normal") W <- rnorm(N, 1, 2)       # W ~ normal covariate.
  if (dist.W == "Binary") W <- rbinom(N, 1, 0.5)    # W ~ binary covariate.
  
  # Relationship between capture probability and covariate.
  
  P_i <- plogis(beta_0 + beta_1*X + beta_2*W)
  
  # Generate a single capture-recapture data set here.
  
  cap.hist0 <- matrix(rbinom(N * tau, 1, P_i), nrow = N, ncol = tau)
  
  cap.hist <- cap.hist0[rowSums(cap.hist0) != 0, ]  # Remove zero (unobserved individuals) rows.
  
  D <- nrow(cap.hist)               # Unique number of individuals captured.
  
  X_i <- X[rowSums(cap.hist0) != 0]   # Covariate X for captured individuals.
  W_i <- W[rowSums(cap.hist0) != 0]   # Covariate W for captured individuals.
  
  X.obs <- cbind(X_i, W_i)
  colnames(X.obs) <- c("x.obs", "w.obs")
  
  # Construct (Binomial) recapture-based quantities under the partial likelihood 
  # for fitting M_h-only models.
  
  pl_data <- prepare_recap_data_Mh(cap.hist, X.obs, tau)
  
  # Fitting M_h-type (GLM) models here.
  
  # Incorrectly specified: GLM with only one covariate X.
  
  model1 <- glm(cbind(Y_recap, n_opp - Y_recap) ~ x.obs,
                      family = binomial(link = "logit"), data = pl_data)
  
  residuals_output1 <- simulateResiduals(model1, n = 1000)  # DHARMa simulated residuals.
  
  # Get population size estimates.
  
  p_i_hat <- predict(model1, newdata = data.frame(x.obs = X_i), type = "response")
  pi_i <- 1 - (1 - p_i_hat)^tau
  N_hat1 <- sum(1/pi_i)     # The Horvitz-Thompson estimator.
  
  # Correctly specified: GLM with an additional covariate.
  
  model2 <- glm(cbind(Y_recap, n_opp - Y_recap) ~ x.obs + w.obs,
                family = binomial(link = "logit"), data = pl_data)
  
  residuals_output2 <- simulateResiduals(model2, n = 1000)
  
  # Get population size estimates.
  
  p_i_hat <- predict(model2, newdata = data.frame(x.obs = X_i, w.obs = W_i), 
                     type = "response")
  pi_i <- 1 - (1 - p_i_hat)^tau
  N_hat2 <- sum(1/pi_i)       # The Horvitz-Thompson estimator.
  
  # Fit model M_h with VGAM (which uses the Huggins' conditional likelihood) 
  # to get the pop size.
  
  data_Mh_VGLM <- data.frame(cbind(cap.hist, X_i, W_i))
  
  if (tau == 7) {
    colnames(data_Mh_VGLM) <- c("y1", "y2", "y3", "y4", "y5", "y6", "y7", 
                                "X_i", "W_i")
  
  M_h_VGLM <- vglm(cbind(y1, y2, y3, y4, y5, y6, y7) ~ X_i + W_i,
                   posbernoulli.t(parallel = TRUE ~ X_i + W_i), 
                   data = data_Mh_VGLM, trace = FALSE)
  }
  
  if (tau == 15) {
    colnames(data_Mh_VGLM) <- c("y1", "y2", "y3", "y4", "y5", 
                                "y6", "y7", "y8", "y9", "y10", 
                                "y11", "y12", "y13", "y14", "y15", "X_i", "W_i")
    
    M_h_VGLM <- vglm(cbind(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, 
                           y11, y12, y13, y14, y15) ~ X_i + W_i,
                     posbernoulli.t(parallel = TRUE ~ X_i + W_i), 
                     data = data_Mh_VGLM, trace = FALSE)
  }
  
  N_hat3 <- M_h_VGLM@extra$N.hat
  
  # DHARMa plots: Side-by-side plots.
  
  figname <- file.path(folder_name, paste0("sim1_GOF_plots_dataset_", ii, ".png"))
  
  png(file = figname, width = 11, height = 9, units = "in", res = 300)
  
  par(mfrow = c(2, 2))
  
  plotQQunif(residuals_output1, 
             main = paste("Dataset_", ii, ": Model 1 (Misspecified)"))
  plotResiduals(residuals_output1)
  plotQQunif(residuals_output2, 
             main = paste("Dataset_", ii, ": Model 2 (Correctly specified)"))
  plotResiduals(residuals_output2)
  
  dev.off()
  
  res[ii, ] <- c(D, N_hat1, N_hat2, N_hat3) 
}

#...............................................................................

res_all <- cbind(apply(res, 2, mean), apply(res, 2, median), 
                 apply(res, 2, sd), apply(res, 2, mad), 
                 sqrt(apply((res - N)^2, 2, mean)))

colnames(res_all) <- c("Mean(N.hat)", "Median(N.hat)", 
                       "SD(N.hat)", "MAD(N.hat)", "RMSE")
rownames(res_all) <- c("D", "PL-GLM [X]", "PL-GLM [X, W]", "CL-GLM [X, W]")
res_all

# Relative % bias for N:

((res_all[, 2] - N)/N)*100
