# ==============================================================================
# CR_GOF_Examples.r
#
# GOF diagnostic using the DHARMa R-package on four real data sets.
#
# A full description of these data and their analysis appears in the 
# Web Sections C1-C4 of the Supplementary Materials
# ==============================================================================

#-------------------------------------------------------------------------------
# Real data example 1: The 2003 Possum data from Stoklosa et al. (2011).
#-------------------------------------------------------------------------------

library(DHARMa)
library(mgcv)
library(tidyverse)
library(secr)

source("CR_GOF_Functions.R")

load("poss2003.RData")

# Remove two outliers as in Huggins and Hwang (2007) and Stoklosa et al. (2011).

body.weight <- scale(data$x.obs[-c(2, 27)])
colnames(body.weight) <- "body.weight"

cap.hist <- data$cap.hist[-c(2, 27), ]

tau <- ncol(cap.hist)    # Number of capture occasions.

D <- nrow(cap.hist)     # Number of unique individuals captured.

#...............................................................................

# Test data for closure using the Stanley--Burnham test via the secr R-package.

captXY <- c()

for(i in 1: nrow(cap.hist)) {
  Occasion <- which (cap.hist[i, ] > 0)
  yy <- length(Occasion)
  ID <- rep(i, yy)
  Session <- rep(1, yy)
  ind1 <- cbind(Session, ID, Occasion)
  captXY <- rbind(captXY, ind1)
}

closure.test(make.capthist(captXY, traps = NULL, 
                           fmt = c("trapID"), noccasions = tau))

# A low p-value (typically less than 0.05) suggests that the population is NOT 
# closed and that an open population model (like the Jolly-Seber, CJS, POPAN 
# model) should be used instead.

# Here the P-value is 0.53320675, it's likely to be a closed population.

#...............................................................................

# Construct (Binomial) recapture-based quantities under the partial likelihood 
# for fitting M_h-only models.

pl_data <- prepare_recap_data_Mh(cap.hist, body.weight, tau)

# Fit a M_h GLM.

mod_mh_glm <- glm(cbind(Y_recap, n_opp - Y_recap) ~ body.weight, 
                 family = binomial(link = "logit"), data = pl_data)

residuals_output1 <- simulateResiduals(mod_mh_glm, n = 2000)

pdf("Example1_GLM.pdf", width = 13, height = 9)
plot(residuals_output1)
dev.off()

# Fit a M_h GAM.

mod_mh_gam <- gam(cbind(Y_recap, n_opp - Y_recap) ~ s(body.weight, bs = "ps"), 
                family = binomial(link = "logit"), gamma = 1.4, data = pl_data)

residuals_output2 <- simulateResiduals(mod_mh_gam, n = 2000)

# Construct (Bernoulli) recapture-based quantities under the partial likelihood
# for fitting M_t- and M_b-type models.

pl_data <- prepare_recap_data_Mtbh(cap.hist, body.weight, tau)

# Fit an M_bh GAM (accounts for behave effects and hete.).

mod_mbh_gam <- gam(Y_recap ~ s(body.weight, bs = "ps") + cap_prev,
               family = binomial(link = "logit"), gamma = 1.4, data = pl_data)

par(mfrow = c(1, 2))
residuals_output3 <- simulateResiduals(mod_mbh_gam, n = 2000)
res_grouped3 <- recalculateResiduals(residuals_output3, group = pl_data$ID)
group_weight <- aggregate(body.weight ~ ID, data = pl_data, FUN = mean)$body.weight
plotQQunif(res_grouped3)
plotResiduals(res_grouped3, group_weight)

# Fit a M_tbh-type GAM (accounts for time, behave effects and hete.).

mod_mtbh_gam <- gam(Y_recap ~ s(body.weight, bs = "ps") + cap_prev + occasion, 
                family = binomial(link = "logit"), gamma = 1.4, data = pl_data)

par(mfrow = c(1, 2))
residuals_output4 <- simulateResiduals(mod_mtbh_gam, n = 2000)
res_grouped4 <- recalculateResiduals(residuals_output4, group = pl_data$ID)
plotQQunif(res_grouped4)
plotResiduals(res_grouped4, group_weight)

model_names <- c("(a) Model: M_h GLM", "(b) Model: M_h GAM", 
                 "(c) Model: M_bh GAM", "(d) Model: M_tbh GAM")
models <- list(mod_mh_glm, mod_mh_gam, mod_mbh_gam, mod_mtbh_gam)

for(i in 1:4) {
  pdf(paste0("Example1_GOF_Model_", i, ".pdf"), width = 14, height = 7)
  if(i == 1 | i == 2) {
    plot(DHARMa::simulateResiduals(models[[i]], n = 2000,
                                          main = paste("Model:", model_names[i])))
  }
  if(i == 3 | i == 4) {
    
    group_weight <- aggregate(body.weight ~ ID, data = pl_data, 
                              FUN = mean)$body.weight
    par(mfrow = c(1, 2))
    res_plot0 <- DHARMa::simulateResiduals(models[[i]], n = 2000,
                                            main = paste("Model:", 
                                                         model_names[i]))
    res_plot <- recalculateResiduals(res_plot0, group = pl_data$ID)
    
    plotQQunif(res_plot)
    plotResiduals(res_plot, group_weight)
  }
  
  mtext(paste("", model_names[i]), side = 3, line = -1.5, 
        outer = TRUE, adj = 0.1, cex = 1.7)
  dev.off()
}

#...............................................................................

# Get the population size from fitting model a M_bh GAM using the partial likelihood
# and the Huggins' conditional likelihood (CL).

mbh_study_grid <- expand.grid(ID = unique(pl_data$ID), occ = 1:tau)

# For Huggins M_bh, we usually assume cap_prev = 0 for the "potential" 
# probability of the first capture.

mbh_study_grid$body.weight <- body.weight[mbh_study_grid$ID] 
mbh_study_grid$cap_prev <- 0 

mbh_preds <- predict(mod_mbh_gam, newdata = mbh_study_grid, type = "response")
mbh_preds_list <- split(mbh_preds, mbh_study_grid$ID)

pi_i <- sapply(mbh_preds_list, function(x) 1 - prod(1 - x))
pi_i2 <- sapply(mbh_preds_list, function(x) prod(1 - x))

N_hat1 <- sum(1/pi_i)         # The Horvitz-Thompson estimator.

# To obtain the standard error for the N_hat1, we use a Taylor's series approx.
# see Stoklosa et al. (2011) for further details.

var.beta <- mod_mbh_gam$Ve
X_des <- model.matrix(mod_mbh_gam, newdata = mbh_study_grid)

pi02 <- (1 - pi_i)/pi_i/pi_i
cz0.dot <- absorb(X_des*c(mbh_preds), tau)
d1 <- apply(pi02*cz0.dot, 2, sum)
var_N_hat1 <- sum((1 - pi_i)/pi_i/pi_i) + d1%*%var.beta%*%d1

suppressWarnings(detach("package:mgcv", unload = TRUE))

library(VGAM)

data_Mbh_VGAM <- data.frame(cbind(cap.hist, body.weight))

colnames(data_Mbh_VGAM) <- c("y1", "y2", "y3", "y4", "y5", "body.weight")
  
M_bh_VGAM <- vgam(cbind(y1, y2, y3, y4, y5) ~ s(body.weight, df = 1.5), 
               posbernoulli.b, data = data_Mbh_VGAM, trace = FALSE) 

N_hat2 <- M_bh_VGAM@extra$N.hat

# Population size estimate (standard error) for M_bh using partial likelihood:

c(N_hat1, sqrt(var_N_hat1))

# Population size estimate (standard error) for M_bh using conditional likelihood:

c(N_hat2, M_bh_VGAM@extra$SE.N.hat) # Standard errors for N_hat.

#-------------------------------------------------------------------------------
# Real data example 2: The Reid deer mice data from Otis et al. (1978).
# We obtained the deermice dataset from the VGAM R-package.
#-------------------------------------------------------------------------------

library(DHARMa)
library(VGAM)
library(secr)

source("CR_GOF_Functions.R")

tau <- 6               # Number of capture occasions.

cap.hist <- deermice[, 1:tau]

D <- NROW(cap.hist)    # Number of unique individuals captured.

age <- as.numeric(deermice$age == "a")   # Extract the age covariate.
sex <- deermice$sex
weight <- deermice$weight

#...............................................................................

# Test data for closure using the Stanley--Burnham test via the secr R-package.

captXY <- c()

for(i in 1: nrow(cap.hist)) {
  Occasion <- which (cap.hist[i, ] > 0)
  yy <- length(Occasion)
  ID <- rep(i, yy)
  Session <- rep(1, yy)
  ind1 <- cbind(Session, ID, Occasion)
  captXY <- rbind(captXY, ind1)
}

closure.test(make.capthist(captXY, traps = NULL, 
                           fmt = c("trapID"), noccasions = tau))

# A low p-value (typically less than 0.05) suggests that the population is NOT 
# closed and that an open population model (like the Jolly-Seber, CJS, POPAN 
# model) should be used instead.

# Here, the P-value is 0.7778398, it's likely to be a closed population.

#...............................................................................

# Construct (Bernoulli) recapture-based quantities under the partial likelihood
# for fitting M_t- and M_b-type models.

X.obs <- cbind(age, sex, weight)
colnames(X.obs) <- c("age", "sex", "weight")

pl_data <- prepare_recap_data_Mtbh(as.matrix(cap.hist), X.obs, tau)

# Fit a M_b GLM (accounts for behavior effects only).

mod_mb <- glm(Y_recap ~ cap_prev, family = binomial(link = "logit"), 
              data = pl_data)

residuals_output1 <- simulateResiduals(mod_mb, n = 1000)
res_grouped1 <- recalculateResiduals(residuals_output1, group = pl_data$ID)

# Fit a M_bh-type GLM (accounts for behavior effects and hete.).

mod_mbh <- glm(Y_recap ~ age + sex + weight + cap_prev, 
                family = binomial(link = "logit"), data = pl_data)

residuals_output2 <- simulateResiduals(mod_mbh, n = 1000)
res_grouped2 <- recalculateResiduals(residuals_output2, group = pl_data$ID)

# Fit a M_tbh-type GLM (accounts for behave effects, time, and hete.).

mod_mtbh <- glm(Y_recap ~ age + sex + weight + cap_prev + occasion, 
               family =  binomial(link = "logit"), data = pl_data)

residuals_output3 <- simulateResiduals(mod_mtbh, n = 1000)
res_grouped3 <- recalculateResiduals(residuals_output3, group = pl_data$ID)

model_names <- c("(a) Model: M_b GLM", "(b) M_bh GLM", "(c) Model: M_tbh GLM")
models <- list(mod_mb, mod_mbh, mod_mtbh)

for(i in 1:3) {
  pdf(paste0("Example2_GOF_Model_", i, ".pdf"), width = 10, height = 5)
  res_plot <- DHARMa::simulateResiduals(models[[i]], n = 1000,
                                        main = paste("Model:", model_names[i]))
  plot(recalculateResiduals(res_plot, group = pl_data$ID))
  
  mtext(paste("", model_names[i]), side = 3, line = -2, 
        outer = TRUE, adj = 0.1, cex = 1.7)
  dev.off()
}

#-------------------------------------------------------------------------------
# Real data example 3: The (well-known) taxicab data from Carothers (1973).
# We obtained the taxicab dataset (Scheme A) from the MARK software White & 
# Burnham (1999). 
#-------------------------------------------------------------------------------

library(RMark)
library(secr)
library(DHARMa)
library(TMB)
library(tidyverse)
library(VGAM)

source("CR_GOF_Functions.R")

# The code below is required for fitting the M_th GLMM. It is a fully marginalized 
# CL approach that uses Laplace approximation via TMB for parameter estimation.

filename_RE <- "TMB_eps_NEW2.cpp"
modelname_RE <- "TMB_eps_NEW2"
TMB::compile(filename_RE, 
             flags = "-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign")

dyn.load(dynlib(modelname_RE))

filepath <- "TaxiM.inp"
lines <- readLines(filepath)
data_lines <- lines[grepl(";", lines)]
ch_strings <- gsub(".*\\*/\\s*([01]+).*", "\\1", data_lines)
ch_list <- strsplit(ch_strings, "")
cap.hist <- matrix(as.numeric(unlist(ch_list)), ncol = nchar(ch_strings[1]), 
                   byrow = TRUE)

tau <- 10               # Number of capture occasions.

D <- nrow(cap.hist)     # Number of unique individuals captured.

#...............................................................................

# Test data for closure using the Stanley--Burnham test via the secr R-package.

captXY <- c()

for(i in 1: nrow(cap.hist)) {
  Occasion <- which (cap.hist[i, ] > 0)
  yy <- length(Occasion)
  ID <- rep(i, yy)
  Session <- rep(1, yy)
  ind1 <- cbind(Session, ID, Occasion)
  captXY <- rbind(captXY, ind1)
}

closure.test(make.capthist(captXY, traps = NULL, fmt = c("trapID"), 
                           noccasions = tau))

# A low p-value (typically less than 0.05) suggests that the population is NOT 
# closed and that an open population model (like the Jolly-Seber, CJS, POPAN 
# model) should be used instead.

# Here, the P-value is 0.4615, it's likely to be a closed population.

#...............................................................................

# Construct (Bernoulli) recapture-based quantities under the partial likelihood
# for fitting M_t- and M_b-type models.

pl_data <- prepare_recap_data_Mtbh(cap.hist, NULL, tau)

# Fit a M_0 GLM (constant model).

mod_m0_glm <- glm(Y_recap ~ 1, data = pl_data, family = binomial(link = "logit"))

residuals_output1 <- simulateResiduals(mod_m0_glm, n = 2000)
res_grouped1 <- recalculateResiduals(residuals_output1, group = pl_data$ID)

# Fit a M_t GLM (accounts for time effects only).

mod_mt_glm <- glm(Y_recap ~ occasion, data = pl_data, 
                  family = binomial(link = "logit"))

residuals_output2 <- simulateResiduals(mod_mt_glm, n = 2000)
res_grouped2 <- recalculateResiduals(residuals_output2, group = pl_data$ID)

# Fit model M_t with VGAM (which uses the Huggins conditional likelihood) 
# to get the population size estimate

VGAM_data_Mt <- data.frame(cap.hist)
colnames(VGAM_data_Mt) <- c("y1", "y2", "y3", "y4", "y5", 
                            "y6", "y7", "y8", "y9", "y10")

M_t <- vglm(cbind(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10) ~ 1, 
            posbernoulli.t, data = VGAM_data_Mt, trace = FALSE) 

N_hat1 <- M_t@extra$N.hat

model_names <- c("(a) Model: M_0 GLM", "(b) M_t GLM")
models <- list(mod_m0_glm, mod_mt_glm)

for(i in 1:2) {
  pdf(paste0("Example3_GOF_Model_", i, ".pdf"), width = 10, height = 5)
  
  res_plot <- DHARMa::simulateResiduals(models[[i]], n = 1000,
                                        main = paste("Model:", model_names[i]))
  plot(recalculateResiduals(res_plot, group = pl_data$ID))
  
  mtext(paste("", model_names[i]), side = 3, line = -2, 
        outer = TRUE, adj = 0.1, cex = 1.7)
  dev.off()
}

#...............................................................................

# Get the population size estimates from fitting model M_th (GLMM) using 
# Huggins conditional likelihood. This model is a fully marginalized CL model
# that uses Laplace approximation via TMB.

tidbits_data <- list(y = as.matrix(cap.hist))

tidbits_parameters <- list(
  time_effects = numeric(ncol(cap.hist)),
  individual_effects = rep(0, D),
  logsigma_individual_effects = log(1) 
)

tidbits_random <- c("individual_effects")

objs <- MakeADFun(data = tidbits_data, 
                  parameters = tidbits_parameters, 
                  random = tidbits_random, 
                  DLL = modelname_RE, silent = TRUE)

fit_RE_marg <- nlminb(start = objs$par, objective = objs$fn, gradient = objs$gr, 
                    control = list(trace = 0, iter.max = 5000, eval.max = 5000))

fit_RE_marg_results <- sdreport(objs, bias.correct = FALSE, 
                                ignore.parm.uncertainty = FALSE, 
                                skip.delta.method = FALSE)

fit_RE_marg_estimates <- as.list(fit_RE_marg_results, what = "Estimate", 
                                 report = TRUE)
fit_RE_marg_all_results <- summary(fit_RE_marg_results, 
                                   select = "report", p.value = TRUE) %>% 
  as.data.frame() %>% rownames_to_column(var = "Parameters")

# Get the population size estimate.

N_hat2 <- fit_RE_marg_all_results %>% dplyr::filter(Parameters == "population_size")

# Population size estimate (standard error) for M_t:

c(N_hat1, M_t@extra$SE.N.hat) 

# Population size estimate (standard error) for M_th:

c(N_hat2$Estimate, N_hat2$`Std. Error`) 

#-------------------------------------------------------------------------------
# Real data example 4: The Dipper data from Lebreton et al. (1992).
# We obtained the dipper dataset from the marked R-package.
#-------------------------------------------------------------------------------

library(marked)
library(DHARMa)
library(secr)

source("CR_GOF_Functions.R")

# Extract the capture history matrix (as a data frame column).

data(dipper)
capture_history <- dipper$ch

# Extract the sex covariate.

sex <- dipper$sex

ch_matrix <- t(sapply(capture_history, function(x) as.numeric(strsplit(x, "")[[1]])))
colnames(ch_matrix) <- paste0("occasion", 1:ncol(ch_matrix))

cap.hist <- ch_matrix 

tau <- ncol(cap.hist)     # Number of capture occasions.
D <- nrow(cap.hist)       # Number of unique individuals captured.

sex <- as.matrix(sex)
colnames(sex) <- c("sex")

#...............................................................................

# Test data for closure using the Stanley--Burnham test via the secr R-package.

captXY <- c()

for(i in 1: nrow(cap.hist)) {
  Occasion <- which (cap.hist[i, ] > 0)
  yy <- length(Occasion)
  ID <- rep(i, yy)
  Session <- rep(1, yy)
  ind1 <- cbind(Session, ID, Occasion)
  captXY <- rbind(captXY, ind1)
}

closure.test(make.capthist(captXY, traps = NULL, 
                           fmt = c("trapID"), noccasions = tau))

# A low p-value (typically less than 0.05) suggests that the population is NOT 
# closed and that an open population model (like the Jolly-Seber, CJS, POPAN 
# model) should be used instead.

# Here the P-value is 5.6315e-38, suggesting we do NOT have a closed population.

#...............................................................................

# Construct (Bernoulli) recapture-based quantities under the partial likelihood.
# for fitting M_t- and M_b-type models.

pl_data <- prepare_recap_data_Mtbh(cap.hist, sex, tau)

# Fit a M_t GLM (accounts for time effects only).

mod_mt_glm <- glm(Y_recap ~ occasion, family = binomial(link = "logit"), 
                  data = pl_data)

residuals_output1 <- simulateResiduals(mod_mt_glm, n = 1000)
res_grouped1 <- recalculateResiduals(residuals_output1, group = pl_data$ID)

# Fit a M_b GLM (accounts for behave effects only).

mod_mb_glm <- glm(Y_recap ~ cap_prev, family = binomial(link = "logit"), 
                  data = pl_data)

residuals_output2 <- simulateResiduals(mod_mb_glm, n = 1000)
res_grouped2 <- recalculateResiduals(residuals_output2, group = pl_data$ID)

# Construct (Binomial) recapture-based quantities under the partial likelihood 
# for fitting M_h models.

pl_data2 <- prepare_recap_data_Mh(cap.hist, sex, tau)

# Fit a M_h GLM.

mod_mh_glm <- glm(cbind(Y_recap, n_opp - Y_recap) ~ sex, 
                  family = binomial(link = "logit"), data = pl_data2)

residuals_output3 <- simulateResiduals(mod_mh_glm, n = 1000)

model_names <- c("(a) Model: M_t GLM", "(b) Model: M_b GLM", "(c) Model: M_h GLM")
models <- list(mod_mt_glm, mod_mb_glm, mod_mh_glm)

for(i in 1:3) {
  pdf(paste0("Example4_GOF_Model_", i, ".pdf"), width = 14, height = 7)
  if(i == 1 | i == 2) {
    res_plot <- DHARMa::simulateResiduals(models[[i]], n = 1000,
                                          main = paste("Model:", model_names[i]))
    plot(recalculateResiduals(res_plot, group = pl_data$ID))
  }
  if(i == 3) plot(DHARMa::simulateResiduals(models[[i]], n = 1000,
                                                     main = paste("Model:", 
                                                                  model_names[i])))
  mtext(paste("", model_names[i]), side = 3, line = -2, 
        outer = TRUE, adj = 0.1, cex = 1.7)
  dev.off()
}
