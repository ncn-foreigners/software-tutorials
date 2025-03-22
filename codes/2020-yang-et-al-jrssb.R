library(nonprobsvy)
library(sampling)
library(doSNOW)
library(progress)
library(foreach)
library(data.table)
library(xtable)

## generate data

seed_for_sim <- 2025
set.seed(seed_for_sim)

N <- 10000
n_A <- 500
sims <- 10 ## ~25 min
p <- 50
KK <- 5
alpha_vec1 <- c(-2, 1, 1, 1, 1, rep(0, p - 5))
alpha_vec2 <- c(0, 0, 0, 3, 3, 3, 3, rep(0, p - 7))
beta_vec <- c(1, 0, 0, 1, 1, 1, 1, rep(0, p - 7))

## generate X
X <- cbind(
  "(Intercept)" = 1, 
  matrix(
    rnorm(N * (p - 1)), 
    nrow = N, byrow = TRUE, 
    dimnames = list(NULL, paste0("X", 1:(p - 1)))
  )
)
X_formula  <- as.formula(paste("~", paste0("X", 1:(p - 1), collapse = " + ")))

## generate Y
Y_11 <- as.numeric(X %*% beta_vec) +  rnorm(N) ## OM I: linear model
Y_12 <- 1 + exp(3*sin(as.numeric(X %*% beta_vec))) + X[, "X5"] + X[, "X6"] + rnorm(N) ## OM II: nonlinear model
pi_Y_21 <- plogis(as.numeric(X %*% beta_vec)) ## OM III: linear model for binary Y
pi_Y_22 <- plogis(2 - log((X %*% beta_vec)^2) - 2*X[,"X5"] + 2*X[, "X6"]) ## OM IV: nonlinear model for binary Y
Y_21 <- rbinom(N, 1, prob = pi_Y_21)
Y_22 <- rbinom(N, 1, prob = pi_Y_22)

## generate probs
pi_A_Y11 <- inclusionprobabilities(0.25 + abs(X[, "X1"]) + 0.03*abs(Y_11), n_A)
pi_A_Y12 <- inclusionprobabilities(0.25 + abs(X[, "X1"]) + 0.03*abs(Y_12), n_A)
pi_A_Y21 <- inclusionprobabilities(0.25 + abs(X[, "X1"]) + 0.03*abs(Y_21), n_A)
pi_A_Y22 <- inclusionprobabilities(0.25 + abs(X[, "X1"]) + 0.03*abs(Y_22), n_A)

pi_B1 <- plogis(as.numeric(X %*% alpha_vec1)) ## PSM I: linear probability
pi_B2 <- plogis(3.5 + as.numeric(log(X^2) %*% alpha_vec2) - sin(X[, "X3"] + X[, "X4"]) - X[,"X5"] + X[, "X6"]) ## PSM II: nonlinear 


## generate data
population <- data.frame(pi_A_Y11, pi_A_Y12, pi_A_Y21, pi_A_Y22, Y_11, Y_12, Y_21, Y_22, X[, 2:p])

controls_inf_paper <- control_inf(vars_selection = TRUE, 
                                  bias_correction = TRUE,
                                  vars_combine = TRUE)

control_sel_paper <- control_sel(nfolds = 5, 
                                 nlambda = 5)

control_sel_paper_h <- control_sel(nfolds = 5, 
                                   nlambda = 5, 
                                   est_method = "gee",
                                   nleqslv_method = "Newton")

A <- Sys.time()

cores <- 8
cl <- makeCluster(cores)
clusterExport(cl, c("N"))

registerDoSNOW(cl)

pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = sims)

opts <- list(progress = \(n) pb$tick())


res <- foreach(k=1:sims, .combine = rbind,
               .packages = c("survey", "nonprobsvy", "sampling"),
               .errorhandling = "stop",
               .options.snow = opts) %dopar% {
                 ## generate samples
                 flag_B1 <- rbinom(N, 1, prob = pi_B1)
                 flag_B2 <- rbinom(N, 1, prob = pi_B2)
                 
                 flag_A_Y11 <- UPpoisson(pik = pi_A_Y11)
                 flag_A_Y12 <- UPpoisson(pik = pi_A_Y12)
                 flag_A_Y21 <- UPpoisson(pik = pi_A_Y21)
                 flag_A_Y22 <- UPpoisson(pik = pi_A_Y22)
                 
                 sample_A_svy_Y11 <- svydesign(ids = ~ 1, probs = ~ pi_A_Y11, pps = "brewer", data = population[flag_A_Y11 == 1, ])
                 sample_A_svy_Y12 <- svydesign(ids = ~ 1, probs = ~ pi_A_Y12, pps = "brewer", data = population[flag_A_Y12 == 1, ])
                 sample_A_svy_Y21 <- svydesign(ids = ~ 1, probs = ~ pi_A_Y21, pps = "brewer", data = population[flag_A_Y21 == 1, ])
                 sample_A_svy_Y22 <- svydesign(ids = ~ 1, probs = ~ pi_A_Y22, pps = "brewer", data = population[flag_A_Y22 == 1, ])
                 
                 sample_B1 <- population[flag_B1 == 1, ]
                 sample_B2 <- population[flag_B2 == 1, ]
                 
                 ### GLM -----------------------------------------------------------------
                 est_b1_glm_y11 <- nonprob(
                   outcome = as.formula(paste("Y_11", X_formula)),
                   data = sample_B1,
                   svydesign = sample_A_svy_Y11,
                   method_outcome = "glm",
                   family_outcome = "gaussian",
                   control_inference = control_inf(vars_selection = T)
                 )
                 est_b1_glm_y12 <- nonprob(
                   outcome = as.formula(paste("Y_12", X_formula)),
                   data = sample_B1,
                   svydesign = sample_A_svy_Y12,
                   method_outcome = "glm",
                   family_outcome = "gaussian",
                   control_inference = control_inf(vars_selection = T)
                 )
                 est_b1_glm_y21 <- nonprob(
                   outcome = as.formula(paste("Y_21", X_formula)),
                   data = sample_B1,
                   svydesign = sample_A_svy_Y21,
                   method_outcome = "glm",
                   family_outcome = "binomial",
                   control_inference = control_inf(vars_selection = T)
                 )
                 est_b1_glm_y22 <- nonprob(
                   outcome = as.formula(paste("Y_22", X_formula)),
                   data = sample_B1,
                   svydesign = sample_A_svy_Y22,
                   method_outcome = "glm",
                   family_outcome = "binomial",
                   control_inference = control_inf(vars_selection = T)
                 )
                 ### 
                 est_b2_glm_y11 <- nonprob(
                   outcome = as.formula(paste("Y_11", X_formula)),
                   data = sample_B2,
                   svydesign = sample_A_svy_Y11,
                   method_outcome = "glm",
                   family_outcome = "gaussian",
                   control_inference = control_inf(vars_selection = T)
                 )
                 est_b2_glm_y12 <- nonprob(
                   outcome = as.formula(paste("Y_12", X_formula)),
                   data = sample_B2,
                   svydesign = sample_A_svy_Y12,
                   method_outcome = "glm",
                   family_outcome = "gaussian",
                   control_inference = control_inf(vars_selection = T)
                 )
                 est_b2_glm_y21 <- nonprob(
                   outcome = as.formula(paste("Y_21", X_formula)),
                   data = sample_B2,
                   svydesign = sample_A_svy_Y21,
                   method_outcome = "glm",
                   family_outcome = "binomial",
                   control_inference = control_inf(vars_selection = T)
                 )
                 est_b2_glm_y22 <- nonprob(
                   outcome = as.formula(paste("Y_22", X_formula)),
                   data = sample_B2,
                   svydesign = sample_A_svy_Y22,
                   method_outcome = "glm",
                   family_outcome = "binomial",
                   control_inference = control_inf(vars_selection = T)
                 )
                              
                 ### IPW -----------------------------------------------------------------
                 est_b1_ipw_y11 <- nonprob(selection = X_formula,
                                           target = ~ Y_11,
                                           data = sample_B1,
                                           svydesign = sample_A_svy_Y11,
                                           control_selection = control_sel_paper,
                                           control_inference = control_inf(vars_selection = TRUE))
                 est_b1_ipw_y12 <- nonprob(selection = X_formula,
                                           target = ~ Y_12,
                                           data = sample_B1,
                                           svydesign = sample_A_svy_Y12,
                                           control_selection = control_sel_paper,
                                           control_inference = control_inf(vars_selection = TRUE))
                 est_b1_ipw_y21 <- nonprob(selection = X_formula,
                                           target = ~ Y_21,
                                           data = sample_B1,
                                           svydesign = sample_A_svy_Y21,
                                           control_selection = control_sel_paper,
                                           control_inference = control_inf(vars_selection = TRUE))
                 est_b1_ipw_y22 <- nonprob(selection = X_formula,
                                           target = ~ Y_22,
                                           data = sample_B1,
                                           svydesign = sample_A_svy_Y22,
                                           control_selection = control_sel_paper,
                                           control_inference = control_inf(vars_selection = TRUE))
                 ### 
                 est_b2_ipw_y11 <- nonprob(selection = X_formula,
                                           target = ~ Y_11,
                                           data = sample_B2,
                                           svydesign = sample_A_svy_Y11,
                                           control_selection = control_sel_paper,
                                           control_inference = control_inf(vars_selection = TRUE))
                 est_b2_ipw_y12 <- nonprob(selection = X_formula,
                                           target = ~ Y_12,
                                           data = sample_B2,
                                           svydesign = sample_A_svy_Y12,
                                           control_selection = control_sel_paper,
                                           control_inference = control_inf(vars_selection = TRUE))
                 est_b2_ipw_y21 <- nonprob(selection = X_formula,
                                           target = ~ Y_21,
                                           data = sample_B2,
                                           svydesign = sample_A_svy_Y21,
                                           control_selection = control_sel_paper,
                                           control_inference = control_inf(vars_selection = TRUE))
                 est_b2_ipw_y22 <- nonprob(selection = X_formula,
                                           target = ~ Y_22,
                                           data = sample_B2,
                                           svydesign = sample_A_svy_Y22,
                                           control_selection = control_sel_paper,
                                           control_inference = control_inf(vars_selection = TRUE))
                 
                 ### DR -----------------------------------------------------------------
                 est_b1_dr_y11 <- nonprob(
                   outcome = as.formula(paste("Y_11", X_formula)),
                   selection = X_formula,
                   data = sample_B1,
                   svydesign = sample_A_svy_Y11,
                   method_outcome = "glm",
                   family_outcome = "gaussian",
                   control_selection = control_sel_paper_h,
                   control_inference = controls_inf_paper,
                   pop_size = N
                 )
                 est_b1_dr_y12 <- nonprob(
                   outcome = as.formula(paste("Y_12", X_formula)),
                   selection = X_formula,
                   data = sample_B1,
                   svydesign = sample_A_svy_Y12,
                   method_outcome = "glm",
                   family_outcome = "gaussian",
                   control_selection = control_sel_paper_h,
                   control_inference = controls_inf_paper,
                   pop_size = N
                 )
                 ## check this one
                 est_b1_dr_y21 <- nonprob(
                   outcome = as.formula(paste("Y_21", X_formula)),
                   selection = X_formula,
                   data = sample_B1,
                   svydesign = sample_A_svy_Y21,
                   method_outcome = "glm",
                   family_outcome = "binomial",
                   control_selection = control_sel_paper_h,
                   control_inference = controls_inf_paper,
                   pop_size = N
                 )
                 est_b1_dr_y22 <- nonprob(
                   outcome = as.formula(paste("Y_22", X_formula)),
                   selection = X_formula,
                   data = sample_B1,
                   svydesign = sample_A_svy_Y22,
                   method_outcome = "glm",
                   family_outcome = "binomial",
                   control_selection = control_sel_paper_h,
                   control_inference = controls_inf_paper,
                   pop_size = N
                 )
                 ### 
                 est_b2_dr_y11 <- nonprob(
                   outcome = as.formula(paste("Y_11", X_formula)),
                   data = sample_B2,
                   svydesign = sample_A_svy_Y11,
                   method_outcome = "glm",
                   family_outcome = "gaussian",
                   control_selection = control_sel_paper_h,
                   control_inference = controls_inf_paper,
                   pop_size = N
                 )
                 est_b2_dr_y12 <- nonprob(
                   outcome = as.formula(paste("Y_12", X_formula)),
                   data = sample_B2,
                   svydesign = sample_A_svy_Y12,
                   method_outcome = "glm",
                   family_outcome = "gaussian",
                   control_selection = control_sel_paper_h,
                   control_inference = controls_inf_paper,
                   pop_size = N
                 )
                 est_b2_dr_y21 <- nonprob(
                   outcome = as.formula(paste("Y_21", X_formula)),
                   data = sample_B2,
                   svydesign = sample_A_svy_Y21,
                   method_outcome = "glm",
                   family_outcome = "binomial",
                   control_selection = control_sel_paper_h,
                   control_inference = controls_inf_paper,
                   pop_size = N
                 )
                 est_b2_dr_y22 <- nonprob(
                   outcome = as.formula(paste("Y_22", X_formula)),
                   data = sample_B2,
                   svydesign = sample_A_svy_Y22,
                   method_outcome = "glm",
                   family_outcome = "binomial",
                   control_selection = control_sel_paper_h,
                   control_inference = controls_inf_paper,
                   pop_size = N
                 )
                 
                 
                 data.frame(
                   k = k,
                   sample = rep(c("b1", "b2"), each = 4),
                   y = rep(c("y11", "y12", "y21", "y22"), times = 2),
                   trues = rep(
                     c(mean(population$Y_11), mean(population$Y_12), 
                       mean(population$Y_21), mean(population$Y_22)), 2
                   ),
                   naive =  c(
                     mean(sample_B1$Y_11),mean(sample_B1$Y_12),
                     mean(sample_B1$Y_21),mean(sample_B1$Y_22),
                     mean(sample_B2$Y_11),mean(sample_B2$Y_12),
                     mean(sample_B2$Y_21),mean(sample_B2$Y_22)
                   ),
                   glm = c(
                     est_b1_glm_y11$output$mean, est_b1_glm_y12$output$mean, 
                     est_b1_glm_y21$output$mean, est_b1_glm_y22$output$mean,
                     est_b2_glm_y11$output$mean, est_b2_glm_y12$output$mean, 
                     est_b2_glm_y21$output$mean, est_b2_glm_y22$output$mean
                   ),
                   ipw = c(
                     est_b1_ipw_y11$output$mean, est_b1_ipw_y12$output$mean, 
                     est_b1_ipw_y21$output$mean, est_b1_ipw_y22$output$mean,
                     est_b2_ipw_y11$output$mean, est_b2_ipw_y12$output$mean, 
                     est_b2_ipw_y21$output$mean, est_b2_ipw_y22$output$mean
                   ),
                   dr = c(
                     est_b1_dr_y11$output$mean, est_b1_dr_y12$output$mean, 
                     est_b1_dr_y21$output$mean, est_b1_dr_y22$output$mean,
                     est_b2_dr_y11$output$mean, est_b2_dr_y12$output$mean, 
                     est_b2_dr_y21$output$mean, est_b2_dr_y22$output$mean
                   ),
                   glm_ci = c(
                     est_b1_glm_y11$confidence_interval[, 1] < mean(population$Y_11) & mean(population$Y_11) < est_b1_glm_y11$confidence_interval[, 2],
                     est_b1_glm_y12$confidence_interval[, 1] < mean(population$Y_12) & mean(population$Y_12) < est_b1_glm_y12$confidence_interval[, 2],
                     est_b1_glm_y21$confidence_interval[, 1] < mean(population$Y_21) & mean(population$Y_21) < est_b1_glm_y21$confidence_interval[, 2],
                     est_b1_glm_y22$confidence_interval[, 1] < mean(population$Y_22) & mean(population$Y_22) < est_b1_glm_y22$confidence_interval[, 2],
                     est_b2_glm_y11$confidence_interval[, 1] < mean(population$Y_11) & mean(population$Y_11) < est_b2_glm_y11$confidence_interval[, 2],
                     est_b2_glm_y12$confidence_interval[, 1] < mean(population$Y_12) & mean(population$Y_12) < est_b2_glm_y12$confidence_interval[, 2],
                     est_b2_glm_y21$confidence_interval[, 1] < mean(population$Y_21) & mean(population$Y_21) < est_b2_glm_y21$confidence_interval[, 2],
                     est_b2_glm_y22$confidence_interval[, 1] < mean(population$Y_22) & mean(population$Y_22) < est_b2_glm_y22$confidence_interval[, 2]
                   ),
                   ipw_ci = c(
                     est_b1_ipw_y11$confidence_interval[, 1] < mean(population$Y_11) & mean(population$Y_11) < est_b1_ipw_y11$confidence_interval[, 2],
                     est_b1_ipw_y12$confidence_interval[, 1] < mean(population$Y_12) & mean(population$Y_12) < est_b1_ipw_y12$confidence_interval[, 2],
                     est_b1_ipw_y21$confidence_interval[, 1] < mean(population$Y_21) & mean(population$Y_21) < est_b1_ipw_y21$confidence_interval[, 2],
                     est_b1_ipw_y22$confidence_interval[, 1] < mean(population$Y_22) & mean(population$Y_22) < est_b1_ipw_y22$confidence_interval[, 2],
                     est_b2_ipw_y11$confidence_interval[, 1] < mean(population$Y_11) & mean(population$Y_11) < est_b2_ipw_y11$confidence_interval[, 2],
                     est_b2_ipw_y12$confidence_interval[, 1] < mean(population$Y_12) & mean(population$Y_12) < est_b2_ipw_y12$confidence_interval[, 2],
                     est_b2_ipw_y21$confidence_interval[, 1] < mean(population$Y_21) & mean(population$Y_21) < est_b2_ipw_y21$confidence_interval[, 2],
                     est_b2_ipw_y22$confidence_interval[, 1] < mean(population$Y_22) & mean(population$Y_22) < est_b2_ipw_y22$confidence_interval[, 2]
                   ),
                   dr_ci = c(
                     est_b1_dr_y11$confidence_interval[, 1] < mean(population$Y_11) & mean(population$Y_11) < est_b1_dr_y11$confidence_interval[, 2],
                     est_b1_dr_y12$confidence_interval[, 1] < mean(population$Y_12) & mean(population$Y_12) < est_b1_dr_y12$confidence_interval[, 2],
                     est_b1_dr_y21$confidence_interval[, 1] < mean(population$Y_21) & mean(population$Y_21) < est_b1_dr_y21$confidence_interval[, 2],
                     est_b1_dr_y22$confidence_interval[, 1] < mean(population$Y_22) & mean(population$Y_22) < est_b1_dr_y22$confidence_interval[, 2],
                     est_b2_dr_y11$confidence_interval[, 1] < mean(population$Y_11) & mean(population$Y_11) < est_b2_dr_y11$confidence_interval[, 2],
                     est_b2_dr_y12$confidence_interval[, 1] < mean(population$Y_12) & mean(population$Y_12) < est_b2_dr_y12$confidence_interval[, 2],
                     est_b2_dr_y21$confidence_interval[, 1] < mean(population$Y_21) & mean(population$Y_21) < est_b2_dr_y21$confidence_interval[, 2],
                     est_b2_dr_y22$confidence_interval[, 1] < mean(population$Y_22) & mean(population$Y_22) < est_b2_dr_y22$confidence_interval[, 2]
                   )
                 )
               }

stopCluster(cl)
Sys.time() - A

## processing results
setDT(res)

results_simulation1_process <- res |> melt(id.vars = 1:4)
results_simulation1_process[, c("est", "ci"):=tstrsplit(variable, "_")]


saveRDS(results_simulation1_process, file = "results/yang2020-500.rds")

## coverage
tab1 <- results_simulation1_process[is.na(ci), .(bias=(mean(value)-mean(trues))*100, 
                                                 se = sd(value)*100, 
                                                 rmse = sqrt((mean(value)-mean(trues))^2 + var(value))*100), 
                                    keyby=.(sample, y, est)] |>
  melt(id.vars = 1:3) |>
  transform(sample=paste(sample, variable, sep = "_")) |>
  transform(variable=NULL) |>
  dcast(... ~ sample, value.var = "value") |>
  {\(x) x[order(y, est)]}() 

tab2 <- results_simulation1_process[!is.na(ci), .(ci = mean(value)*100), 
                                    keyby=.(sample, y, est)]  |>
  dcast(... ~ sample, value.var = "ci")

tab1 <- tab1[tab2, on = c("y", "est")] |>
  setcolorder(c("y", "est",  "b1_bias", "b1_se", "b1_rmse", "b1", "b2_bias", "b2_se", "b2_rmse", "b2")) 
