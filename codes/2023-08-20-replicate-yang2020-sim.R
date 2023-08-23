# coverage 
library(nonprobsvy)
library(survey)
library(sampling)
library(glue)
library(data.table)
library(ggplot2)
seed_for_sim <- 2023-8-15
set.seed(seed_for_sim)
N <- 10000
n_A <- 500
p <- 50
alpha_vec1 <- c(-2, 1, 1, 1,1, rep(0, p-5))
alpha_vec2 <- c(0,0,0,3,3,3,3, rep(0, p-7))
beta_vec <- c(1,0,0,1,1,1,1, rep(0, p-7))
X <- cbind("(Intercept)"=1, matrix(rnorm(N*(p-1)), nrow=N, byrow=T, dimnames = list(NULL, paste0("X",1:(p-1)))))
X_formula  <- as.formula(paste("~", paste0("X",1:(p-1), collapse = " + ")))
X_totals <- colSums(X)
X_means <- colMeans(X[,-1])
Y_11 <- 1 + as.numeric(X %*% beta_vec) +   rnorm(N) ## OM I: linear model
Y_12 <- 1 + exp(3*sin(as.numeric(X %*% beta_vec))) + X[, "X5"] + X[, "X6"] + rnorm(N) ## OM II: nonlinear model
pi_Y_21 <- plogis(as.numeric(X %*% beta_vec)) ## OM III: linear model for binary Y
pi_Y_22 <- plogis(2 - log((X %*% beta_vec)^2) - 2*X[,"X5"] + 2*X[, "X6"]) ## OM IV: nonlinear model for binary Y
Y_21 <- rbinom(N, 1, prob = pi_Y_21)
Y_22 <- rbinom(N, 1, prob = pi_Y_22)
pi_A <- inclusionprobabilities(0.25 + abs(X[, "X1"]) + 0.03*abs(Y_11), n_A) ## inclusion based on Y_11 only 
pi_B1 <- plogis(as.numeric(X %*% alpha_vec1)) ## PSM I: linear probability
pi_B2 <- plogis(3.5 + as.numeric(log(X^2) %*% alpha_vec2) - sin(X[, "X3"] + X[, "X4"]) - X[,"X5"] + X[, "X6"]) ## PSM II: nonlinear 


results_mi <- list()
## here is sampling
for (b in 1:10) {
  print(b)
  set.seed(b)

  # generate data -----------------------------------------------------------

  flag_B1 <- rbinom(N, 1, prob = pi_B1)
  flag_B2 <- rbinom(N, 1, prob = pi_B2)
  flag_A <- UPpoisson(pik = pi_A)
  pop_data <- data.frame(pi_A, flag_A, flag_B1, flag_B2, Y_11, Y_12, Y_21, Y_22, X[, 2:p])
  sample_A_svy <- svydesign(ids = ~ 1, probs = ~ pi_A, 
                            pps = poisson_sampling(pop_data$pi_A[pop_data$flag_A == 1]), 
                            data = pop_data[pop_data$flag_A == 1, ])
  sample_A_svy_cal <- calibrate(sample_A_svy, 
                                formula = X_formula,
                                population = X_totals, 
                                calfun = cal.raking)
  sample_B1 <- pop_data[pop_data$flag_B1 == 1, ]
  sample_B2 <- pop_data[pop_data$flag_B2 == 1, ]
  
  print("Mass imputation")
  # mass imputation with scad -----------------------------------------------
  mass_imp_glm <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B1) {
    if (x %in% c("Y_11", "Y_12")) {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ ", as.character(X_formula)[2])),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          control_inference = controlInf(vars_selection = TRUE),
                          family_outcome = "gaussian")
    } else {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ ", as.character(X_formula)[2])),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          control_inference = controlInf(vars_selection = TRUE),
                          family_outcome = "binomial")
    }
    data.frame(y=x, est_mi_y$output, est_mi_y$confidence_interval)
   }
  )
  
  mass_imp_glm2 <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B2) {
    if (x %in% c("Y_11", "Y_12")) {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ ", as.character(X_formula)[2])),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          control_inference = controlInf(vars_selection = TRUE),
                          family_outcome = "gaussian")
    } else {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ ", as.character(X_formula)[2])),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          control_inference = controlInf(vars_selection = TRUE),
                          family_outcome = "binomial")
    }
    data.frame(y=x, est_mi_y$output, est_mi_y$confidence_interval)
  }
  )
  
  # mass imputation (nn) with scad -----------------------------------------------
  mass_imp_glm_nn <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), 
                            FUN = function(x, data=sample_B1) {
    if (x %in% c("Y_11", "Y_12")) {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ ", as.character(X_formula)[2])),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          control_inference = controlInf(vars_selection = TRUE),
                          method_outcome = "nn")
    } else {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ ", as.character(X_formula)[2])),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          control_inference = controlInf(vars_selection = TRUE),
                          method_outcome = "nn")
    }
    data.frame(y=x, est_mi_y$output, est_mi_y$confidence_interval)
  }
  )
  
  mass_imp_glm_nn2 <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), 
                             FUN = function(x, data=sample_B2) {
    if (x %in% c("Y_11", "Y_12")) {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ ", as.character(X_formula)[2])),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          control_inference = controlInf(vars_selection = TRUE),
                          family_outcome = "gaussian")
    } else {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ ", as.character(X_formula)[2])),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          control_inference = controlInf(vars_selection = TRUE),
                          family_outcome = "binomial")
    }
    data.frame(y=x, est_mi_y$output, est_mi_y$confidence_interval)
  }
  )
  

  # ipw ---------------------------------------------------------------------
  est_ipw_h1 <- nonprob(selection = X_formula,
                     target = ~ Y_11,
                     data = sample_B1,
                     svydesign = sample_A_svy_cal,
                     control_selection = controlSel(h_x = "1", est_method_sel = "gee"),
                     control_inference = controlInf(vars_selection = TRUE))
  
  est_ipw_h2 <- nonprob(selection = X_formula,
                     target = ~ Y_11,
                     data = sample_B1,
                     svydesign = sample_A_svy_cal,
                     control_selection = controlSel(h_x = "2", est_method_sel = "gee"),
                     control_inference = controlInf(vars_selection = TRUE))

  est_ipw_mle <- nonprob(selection = X_formula,
                        target = ~ Y_11,
                        data = sample_B1,
                        svydesign = sample_A_svy_cal,
                        control_inference = controlInf(vars_selection = TRUE))
  
  # save results ------------------------------------------------------------

  
  ## save results
  results_mi[[b]] <- rbind(
    rbindlist(mass_imp_glm)[, ":="(dataset="PSM I", est="glm")],
    rbindlist(mass_imp_glm2)[, ":="(dataset="PSM II", est="glm")],
    rbindlist(mass_imp_glm_nn)[, ":="(dataset="PSM I", est="nn")],
    rbindlist(mass_imp_glm_nn2)[, ":="(dataset="PSM II", est="nn")],
    data.frame(est_ipw$output,est_ipw$confidence_interval, dataset= "PSM I", est = "ipw (h=1)"),
    data.frame(est_ipw$output,est_ipw$confidence_interval, dataset= "PSM I", est = "ipw (h=2)")
  )
}


results_mi_df <- rbindlist(results_mi, idcol = "b")
results_mi_df[y == "Y_11", true := mean(Y_11)]
results_mi_df[y == "Y_12", true := mean(Y_12)]
results_mi_df[y == "Y_21", true := mean(Y_21)]
results_mi_df[y == "Y_22", true := mean(Y_22)]

## coverage
results_mi_df[, .(m = mean(lower_bound < true & upper_bound > true)), keyby=.(y, dataset, est)]

ggplot(data = results_mi_df, aes(x = dataset, y = mean, fill = est)) + 
  geom_boxplot(position = "dodge") +
  geom_hline(aes(yintercept = true), color = "red", linetype = "dashed") +
  facet_wrap(~ y, nrow = 1, scales = "free_y") 
