# ==============================================================================
# CR_GOF_Simulation_Scenario2.r
#
# Simulation scenario 2: Models where heterogeneity is modeled as a non-linear 
# (smooth) relationship with a single covariate.
#
# We evaluate model GOFs (for 2 models only) using DHARMa. GOF residual plots 
# are assessed visually and we use the KS, Dispersion and Outlier tests. 
# 
# We repeat this 50 times. 
#
# We also report some summaries for population size estimates for each model 
# using partial likelihood with a Horvitz-Thompson estimator. We also include 
# population size estimates fitting a (Huggins) conditional 
# likelihood/Horvitz-Thompson estimator M_h (GAM) using VGAM.
#
# A full description of these data and their analysis appears in the 
# Supplementary Materials (Web Section C2) of the main text.
# ==============================================================================

library(DHARMa)

source("CR_GOF_Functions.R")

nsim <- 50   # Number of simulations.

# Set parameters.

N <- 50       # True population size.
tau <- 7       # Number of capture occasions.

if(tau != 7 && tau != 15) stop("Unsupported number of capture occasions (tau). 
                               Tau needs to be either 7 or 15 for these sims.")

beta_0 <- -1.2
beta_1 <- 2

res <- matrix(NA, nsim, 4)

folder_name <- "Sim2_GOF"
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}

set.seed(123)

for (ii in 1:nsim) {

  X <- rnorm(N, 0, 1)   # Covariate included in the model.
  
  # Non-linear relationship between capture probability and covariate.
  
  P_i <- plogis(beta_0 + beta_1*sin(1.5*X))
  plot(X, P_i)
  
  # Generate a single capture-recapture data set here.
  
  cap.hist0 <- matrix(rbinom(N * tau, 1, P_i), nrow = N, ncol = tau)
  
  cap.hist <- cap.hist0[rowSums(cap.hist0) != 0, ]  # Remove zero (unobserved individuals) rows.
  
  D <- nrow(cap.hist)          # Unique number of individuals captured.
  
  X_i <- X[rowSums(cap.hist0) != 0]    # Covariate for captured individuals.
  
  x.obs <- as.matrix(X_i)
  colnames(x.obs) <- c("x.obs")
  
  # Construct (Binomial) recapture-based quantities under the partial likelihood 
  # for fitting M_h-only models.
  
  pl_data <- prepare_recap_data_Mh(cap.hist, x.obs, tau)
  
  # Start fitting M_h models here.
  
  # # Incorrectly specified: Linear GLM with one covariate.
  
  model1 <- glm(cbind(Y_recap, n_opp - Y_recap) ~ x.obs,
                family = binomial(link = "logit"), data = pl_data)
  
  residuals_output1 <- simulateResiduals(model1, n = 1000)
  
  # Get population size estimates.
  
  p_i_hat <- predict(model1, newdata = data.frame(x.obs = X_i), type = "response")
  pi_i <- 1 - (1 - p_i_hat)^tau
  N_hat1 <- sum(1/pi_i)       # The Horvitz-Thompson estimator.
  
  # Correctly specified: GAM with a smooth term.

  library(mgcv)
  
  model2 <- gam(cbind(Y_recap, n_opp - Y_recap) ~ s(x.obs, bs = "tp"), 
                family = binomial(link = "logit"), data = pl_data)
  
  residuals_output2 <- simulateResiduals(model2, n = 1000)
  
  # Get population size estimates.
  
  p_i_hat <- predict(model2, newdata = data.frame(x.obs = X_i), type = "response")
  pi_i <- 1 - (1 - p_i_hat)^tau
  N_hat2 <- sum(1/pi_i)      
  
  # Fit model M_h with VGAM (which uses the Huggins' conditional likelihood) 
  # to get the population size estimates.
  
  suppressWarnings(detach("package:mgcv", unload = TRUE))
  
  library(VGAM)
  
  data_Mh_VGAM <- data.frame(cbind(cap.hist, X_i))
  
  if (tau == 7) {
    colnames(data_Mh_VGAM) <- c("y1", "y2", "y3", "y4", "y5", "y6", "y7", "X_i")
    
    # Set df = 5 for N = 50, df = 8 for N = 400 for an optimal smoothing parameter
    
    M_h_VGAM <- vgam(cbind(y1, y2, y3, y4, y5, y6, y7) ~ s(X_i, df = 8),
                   posbernoulli.t(parallel = TRUE ~ s(X_i, df = 8)), 
                   data = data_Mh_VGAM, trace = FALSE)
  }
  
  if (tau == 15) {
    colnames(data_Mh_VGAM) <- c("y1", "y2", "y3", "y4", "y5", "y6", "y7", "y8",
                                "y9", "y10", "y11", "y12", "y13", "y14", "y15", "X_i")
    
    M_h_VGAM <- vgam(cbind(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, 
                           y13, y14, y15) ~ s(X_i, df = 7),
                     posbernoulli.t(parallel = TRUE ~ s(X_i, df = 7)), 
                     data = data_Mh_VGAM, trace = FALSE)
  }
  
  N_hat3 <- M_h_VGAM@extra$N.hat
  
  detach("package:VGAM", unload = TRUE)
  
  # DHARMa plots: Side-by-side plots.
  
  figname <- file.path(folder_name, paste0("sim2_GOF_plots_dataset_", ii, ".png"))
  
  png(file = figname, width = 11, height = 9, units = "in", res = 300)
  
  par(mfrow = c(2, 2))
  
  plotQQunif(residuals_output1, 
             main = paste("Dataset_", ii, ": Model 1 (Misspecified, GLM)"))
  plotResiduals(residuals_output1)
  plotQQunif(residuals_output2, 
             main = paste("Dataset_", ii, ": Model 2 (Correctly specified, GAM)"))
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
rownames(res_all) <- c("D", "PL-GLM", "PL-GAM", "CL-VGAM")
res_all

# Relative % bias for N:

((res_all[, 2] - N)/N)*100

