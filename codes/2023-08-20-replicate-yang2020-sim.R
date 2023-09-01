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
Y_11 <- as.numeric(X %*% beta_vec) +  rnorm(N) ## OM I: linear model
Y_12 <- 1 + exp(3*sin(as.numeric(X %*% beta_vec))) + X[, "X5"] + X[, "X6"] + rnorm(N) ## OM II: nonlinear model
pi_Y_21 <- plogis(as.numeric(X %*% beta_vec)) ## OM III: linear model for binary Y
pi_Y_22 <- plogis(2 - log((X %*% beta_vec)^2) - 2*X[,"X5"] + 2*X[, "X6"]) ## OM IV: nonlinear model for binary Y
Y_21 <- rbinom(N, 1, prob = pi_Y_21)
Y_22 <- rbinom(N, 1, prob = pi_Y_22)
pi_A <- inclusionprobabilities(0.25 + abs(X[, "X1"]) + 0.03*abs(Y_11), n_A) ## inclusion based on Y_11 only 
pi_B1 <- plogis(as.numeric(X %*% alpha_vec1)) ## PSM I: linear probability
pi_B2 <- plogis(3.5 + as.numeric(log(X^2) %*% alpha_vec2) - sin(X[, "X3"] + X[, "X4"]) - X[,"X5"] + X[, "X6"]) ## PSM II: nonlinear 

## checking coverage
result_svy <- list()

for (b in 1:500) {
  set.seed(b)
  print(b)
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
  res1 <- svymean(~Y_11+Y_12+Y_21+Y_22, sample_A_svy)
  res1_ci <- as.data.frame(confint(res1))
  res2 <- svymean(~Y_11+Y_12+Y_21+Y_22, sample_A_svy_cal)
  res2_ci <- as.data.frame(confint(res2))
  result_svy[[b]] <- data.frame(
    y = names(res1),
    est1 = as.numeric(res1),
    est1_lb = res1_ci$`2.5 %`,
    est1_ub = res1_ci$`97.5 %`,
    est2 = as.numeric(res2),
    est2_lb = res1_ci$`2.5 %`,
    est2_ub = res1_ci$`97.5 %`
  )
}
result_svy_df <- rbindlist(result_svy, idcol = 'b')
result_svy_df[y == "Y_11", true := mean(Y_11)]
result_svy_df[y == "Y_12", true := mean(Y_12)]
result_svy_df[y == "Y_21", true := mean(Y_21)]
result_svy_df[y == "Y_22", true := mean(Y_22)]
result_svy_df[, .(est1 = mean(est1_lb < true & est1_ub > true)), keyby=.(y)]
result_svy_df[, .(est2 = mean(est2_lb < true & est2_ub > true)), keyby=.(y)]

results_correct_x <- list()
# no selection of X (correct vars) ----------------------------------------
for (b in 1:500) {
  set.seed(b)
  print(b)
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
  # mass imputation ---------------------------------------------------------
  mass_imp_glm <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B1) {
    if (x %in% c("Y_11", "Y_12")) {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          family_outcome = "gaussian")
    } else {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          family_outcome = "binomial")
    }
    data.frame(y=x, est_mi_y$output, est_mi_y$confidence_interval)
  })
  mass_imp_glm_nn <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B1) {
    if (x %in% c("Y_11", "Y_12")) {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          method_outcome = "nn")
    } else {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          method_outcome = "nn")
    }
    data.frame(y=x, est_mi_y$output, est_mi_y$confidence_interval)
  })
  mass_imp_glm2 <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B2) {
    if (x %in% c("Y_11", "Y_12")) {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          family_outcome = "gaussian")
    } else {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          family_outcome = "binomial")
    }
    data.frame(y=x, est_mi_y$output, est_mi_y$confidence_interval)
  })
  mass_imp_glm_nn2 <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B2) {
    if (x %in% c("Y_11", "Y_12")) {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          method_outcome = "nn")
    } else {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          method_outcome = "nn")
    }
    data.frame(y=x, est_mi_y$output, est_mi_y$confidence_interval)
  })
  # ipw B1 ---------------------------------------------------------------------
  est_ipw_h1_0 <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                          target = ~ Y_11,
                          data = sample_B1,
                          svydesign = sample_A_svy_cal,
                          control_selection = controlSel(h = 1, est_method_sel = "gee"))
  est_ipw_h1_1 <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                          target = ~ Y_12,
                          data = sample_B1,
                          svydesign = sample_A_svy_cal,
                          control_selection = controlSel(h = 1, est_method_sel = "gee"))
  est_ipw_h1_2 <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                          target = ~ Y_21,
                          data = sample_B1,
                          svydesign = sample_A_svy_cal,
                          control_selection = controlSel(h= 1, est_method_sel = "gee"))
  est_ipw_h1_3 <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                          target = ~ Y_22,
                          data = sample_B1,
                          svydesign = sample_A_svy_cal,
                          control_selection = controlSel(h=1, est_method_sel = "gee"))
  
  est_ipw_mle <- nonprob(selection =  ~ X1 + X2 + X3 + X4,
                         target = ~ Y_11 + Y_12 + Y_21 + Y_22,
                         data = sample_B1,
                         svydesign = sample_A_svy_cal)
  # ipw B2 ---------------------------------------------------------------------
  est_ipw_h1_0_2 <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                            target = ~ Y_11,
                            data = sample_B2,
                            svydesign = sample_A_svy_cal,
                            control_selection = controlSel(h=1, est_method_sel = "gee"))
  
  est_ipw_h1_1_2 <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                          target = ~ Y_12,
                          data = sample_B2,
                          svydesign = sample_A_svy_cal,
                          control_selection = controlSel(h=1, est_method_sel = "gee"))
  
  est_ipw_h1_2_2 <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                          target = ~ Y_21,
                          data = sample_B2,
                          svydesign = sample_A_svy_cal,
                          control_selection = controlSel(h=1, est_method_sel = "gee"))
  
  est_ipw_h1_3_2 <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                          target = ~ Y_22,
                          data = sample_B2,
                          svydesign = sample_A_svy_cal,
                          control_selection = controlSel(h=1, est_method_sel = "gee"))
  
  est_ipw_mle_2 <- nonprob(selection =  ~ X1 + X2 + X3 + X4,
                         target = ~ Y_11 + Y_12 + Y_21 + Y_22,
                         data = sample_B2,
                         svydesign = sample_A_svy_cal)
  # dr ----------------------------------------------------------------------
  dr_y11 <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B1) {
    if (x %in% c("Y_11", "Y_12")) {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ {X_formula}")),
                    selection = ~ X1 + X2 + X3 + X4,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    family_outcome = "gaussian")
    } else {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                    selection = ~ X1 + X2 + X3 + X4,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    family_outcome = "binomial")
    }
    data.frame(y=x, dr$output, dr$confidence_interval)
  })
  dr_y11_h1 <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B1) {
    if (x %in% c("Y_11", "Y_12")) {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                    selection = ~ X1 + X2 + X3 + X4,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_selection = controlSel(est_method_sel = "gee", h = 1),
                    family_outcome = "gaussian")
    } else {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                    selection = ~ X1 + X2 + X3 + X4,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_selection = controlSel(est_method_sel = "gee", h = 1),
                    family_outcome = "binomial")
    }
    data.frame(y=x, dr$output, dr$confidence_interval)
  })
  dr_y11_2 <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B2) {
    if (x %in% c("Y_11", "Y_12")) {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                    selection = ~ X1 + X2 + X3 + X4,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    family_outcome = "gaussian")
    } else {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                    selection = ~ X1 + X2 + X3 + X4,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    family_outcome = "binomial")
    }
    data.frame(y=x, dr$output, dr$confidence_interval)
  })
  dr_y11_h1_2 <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B2) {
    if (x %in% c("Y_11", "Y_12")) {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                    selection = ~ X1 + X2 + X3 + X4,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_selection = controlSel(est_method_sel = "gee", h= 1),
                    family_outcome = "gaussian")
    } else {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                    selection = ~ X1 + X2 + X3 + X4,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_selection = controlSel(est_method_sel = "gee", h=1),
                    family_outcome = "binomial")
    }
    data.frame(y=x, dr$output, dr$confidence_interval)
  })
  
  # dr (MM) ----------------------------------------------------------------------
  dr_y11_mm <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B1) {
    if (x %in% c("Y_11", "Y_12")) {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                    selection = ~ X1 + X2 + X3 + X4,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    family_outcome = "gaussian",
                    control_inference = controlInf(bias_correction = TRUE))
    } else {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                    selection = ~ X1 + X2 + X3 + X4,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    family_outcome = "binomial",
                    control_inference = controlInf(bias_correction = TRUE))
    }
    data.frame(y=x, dr$output, dr$confidence_interval)
  })
  dr_y11_h1_mm <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B1) {
    if (x %in% c("Y_11", "Y_12")) {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                    selection = ~ X1 + X2 + X3 + X4,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_selection = controlSel(h = 1),
                    family_outcome = "gaussian",
                    control_inference = controlInf(bias_correction = TRUE))
    } else {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                    selection = ~ X1 + X2 + X3 + X4,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_selection = controlSel(h=1),
                    family_outcome = "binomial",
                    control_inference = controlInf(bias_correction = TRUE))
    }
    data.frame(y=x, dr$output, dr$confidence_interval)
  })
  dr_y11_mm_2 <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B2) {
    if (x %in% c("Y_11", "Y_12")) {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                    selection = ~ X1 + X2 + X3 + X4,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_selection = controlSel(h=2),
                    family_outcome = "gaussian",
                    control_inference = controlInf(bias_correction = TRUE))
    } else {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                    selection = ~ X1 + X2 + X3 + X4,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_selection = controlSel(h=2),
                    family_outcome = "binomial",
                    control_inference = controlInf(bias_correction = TRUE))
    }
    data.frame(y=x, dr$output, dr$confidence_interval)
  })
  dr_y11_h1_mm_2 <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B2) {
    if (x %in% c("Y_11", "Y_12")) {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                    selection = ~ X1 + X2 + X3 + X4,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_selection = controlSel(h=1),
                    family_outcome = "gaussian",
                    control_inference = controlInf(bias_correction = TRUE))
    } else {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                    selection = ~ X1 + X2 + X3 + X4,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_selection = controlSel(h=1),
                    family_outcome = "binomial",
                    control_inference = controlInf(bias_correction = TRUE))
    }
    data.frame(y=x, dr$output, dr$confidence_interval)
  })
  # save results ------------------------------------------------------------
  results_correct_x[[b]] <- rbind(
    rbindlist(mass_imp_glm)[, ":="(dataset="PSM I", est="glm")],
    rbindlist(mass_imp_glm2)[, ":="(dataset="PSM II", est="glm")],
    rbindlist(mass_imp_glm_nn)[, ":="(dataset="PSM I", est="nn")],
    rbindlist(mass_imp_glm_nn2)[, ":="(dataset="PSM II", est="nn")],
    ## PSM I
    data.frame(y=rownames(est_ipw_h1_0$output),
               est_ipw_h1_0$output,est_ipw_h1_0$confidence_interval, dataset= "PSM I", est = "ipw (h=1)"),
    data.frame(y=rownames(est_ipw_h1_1$output),
               est_ipw_h1_1$output,est_ipw_h1_1$confidence_interval, dataset= "PSM I", est = "ipw (h=1)"),
    data.frame(y=rownames(est_ipw_h1_2$output),
               est_ipw_h1_2$output,est_ipw_h1_2$confidence_interval, dataset= "PSM I", est = "ipw (h=1)"),
    data.frame(y=rownames(est_ipw_h1_3$output),
               est_ipw_h1_3$output,est_ipw_h1_3$confidence_interval, dataset= "PSM I", est = "ipw (h=1)"),
    data.frame(y=rownames(est_ipw_mle$output),
               est_ipw_mle$output,est_ipw_mle$confidence_interval, dataset= "PSM I", est = "ipw (h=2)"),
    ## PSM II
    data.frame(y=rownames(est_ipw_h1_0_2$output),
               est_ipw_h1_0_2$output,est_ipw_h1_0_2$confidence_interval, dataset= "PSM II", est = "ipw (h=1)"),
    data.frame(y=rownames(est_ipw_h1_1_2$output),
               est_ipw_h1_1_2$output,est_ipw_h1_1_2$confidence_interval, dataset= "PSM II", est = "ipw (h=1)"),
    data.frame(y=rownames(est_ipw_h1_2_2$output),
               est_ipw_h1_2_2$output,est_ipw_h1_2_2$confidence_interval, dataset= "PSM II", est = "ipw (h=1)"),
    data.frame(y=rownames(est_ipw_h1_3_2$output),
               est_ipw_h1_3_2$output,est_ipw_h1_3_2$confidence_interval, dataset= "PSM II", est = "ipw (h=1)"),
    data.frame(y=rownames(est_ipw_mle_2$output),
               est_ipw_mle_2$output,est_ipw_mle_2$confidence_interval, dataset= "PSM II", est = "ipw (h=2)"),
    ## dr
    rbindlist(dr_y11)[, ":="(dataset="PSM I", est="dr (h=2)")],
    rbindlist(dr_y11_h1)[, ":="(dataset="PSM I", est="dr (h=1)")],
    rbindlist(dr_y11_2)[, ":="(dataset="PSM II", est="dr (h=2)")],
    rbindlist(dr_y11_h1_2)[, ":="(dataset="PSM II", est="dr (h=1)")],
    ## dr (mm)
    rbindlist(dr_y11_mm)[, ":="(dataset="PSM I", est="dr mm (h=2)")],
    rbindlist(dr_y11_h1_mm)[, ":="(dataset="PSM I", est="dr mm (h=1)")],
    rbindlist(dr_y11_mm_2)[, ":="(dataset="PSM II", est="dr mm (h=2)")],
    rbindlist(dr_y11_h1_mm_2)[, ":="(dataset="PSM II", est="dr mm (h=1)")]
  )
}

## collect errors
results_correct_x_df <- rbindlist(results_correct_x, idcol = "b")
results_correct_x_df[y == "Y_11", true := mean(Y_11)]
results_correct_x_df[y == "Y_12", true := mean(Y_12)]
results_correct_x_df[y == "Y_21", true := mean(Y_21)]
results_correct_x_df[y == "Y_22", true := mean(Y_22)]
saveRDS(results_correct_x_df, file = "results/yang2020-correct-x-500.rds")

## coverage
results_correct_x_df[!is.na(lower_bound), 
              .(m = mean(lower_bound < true & upper_bound > true)), keyby=.(y, dataset, est)]

ggplot(data = results_correct_x_df[!is.na(lower_bound) & lower_bound > 0], 
       aes(x = est, y = mean)) + 
  geom_boxplot(position = "dodge") +
  geom_hline(aes(yintercept = true), color = "red", linetype = "dashed") +
  facet_wrap(dataset~y, ncol = 4, scales = "free_y") +
  theme(axis.text.x=element_text(angle=45, vjust=1,hjust=1))


results_correct_x_df[y == "Y_11"
   , .(bias = (mean(mean) - mean(true))/mean(true)*100,
       var = var(mean),
       rmse = sqrt( (mean(mean) - mean(true))^2 + var(mean))), 
   keyby=.(y, dataset, est)]


# selection of variables --------------------------------------------------

X_formula  <- as.formula(paste("~", paste0("X",1:10, collapse = " + ")))
X_totals <- X_totals[1:11]

results_var_sel <- list()
# no selection of X (correct vars) ----------------------------------------
for (b in 1:1) {
  set.seed(b)
  print(b)
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
  # mass imputation ---------------------------------------------------------
  mass_imp_glm <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B1) {
    if (x %in% c("Y_11", "Y_12")) {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          family_outcome = "gaussian",
                          control_inference = controlInf(vars_selection = TRUE))
    } else {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          family_outcome = "binomial",
                          control_inference = controlInf(vars_selection = TRUE))
    }
    data.frame(y=x, est_mi_y$output, est_mi_y$confidence_interval)
  })
  mass_imp_glm_nn <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B1) {
    if (x %in% c("Y_11", "Y_12")) {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          method_outcome = "nn",
                          control_inference = controlInf(vars_selection = TRUE))
    } else {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          method_outcome = "nn",
                          control_inference = controlInf(vars_selection = TRUE))
    }
    data.frame(y=x, est_mi_y$output, est_mi_y$confidence_interval)
  })
  mass_imp_glm2 <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B2) {
    if (x %in% c("Y_11", "Y_12")) {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          family_outcome = "gaussian",
                          control_inference = controlInf(vars_selection = TRUE))
    } else {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          family_outcome = "binomial",
                          control_inference = controlInf(vars_selection = TRUE))
    }
    data.frame(y=x, est_mi_y$output, est_mi_y$confidence_interval)
  })
  mass_imp_glm_nn2 <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B2) {
    if (x %in% c("Y_11", "Y_12")) {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          method_outcome = "nn",
                          control_inference = controlInf(vars_selection = TRUE))
    } else {
      est_mi_y <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                          data = data,
                          svydesign = sample_A_svy_cal,
                          method_outcome = "nn",
                          control_inference = controlInf(vars_selection = TRUE))
    }
    data.frame(y=x, est_mi_y$output, est_mi_y$confidence_interval)
  })
  # ipw B1 ---------------------------------------------------------------------
  est_ipw_h1_0 <- nonprob(selection = X_formula,
                          target = ~ Y_11,
                          data = sample_B1,
                          svydesign = sample_A_svy_cal,
                          verbose = TRUE,
                          control_selection = controlSel(h = 1, est_method_sel = "gee", nfolds = 5),
                          control_inference = controlInf(vars_selection = TRUE))
  
  est_ipw_h1_1 <- nonprob(selection = X_formula,
                          target = ~ Y_12,
                          data = sample_B1,
                          verbose = TRUE,
                          svydesign = sample_A_svy_cal,
                          control_selection = controlSel(h = 1, est_method_sel = "gee", nfolds = 5),
                          control_inference = controlInf(vars_selection = TRUE))
  
  est_ipw_h1_2 <- nonprob(selection = X_formula,
                          target = ~ Y_21,
                          data = sample_B1,
                          verbose = TRUE,
                          svydesign = sample_A_svy_cal,
                          control_selection = controlSel(h= 1, est_method_sel = "gee", nfolds = 5),
                          control_inference = controlInf(vars_selection = TRUE))
  
  est_ipw_h1_3 <- nonprob(selection = X_formula,
                          target = ~ Y_22,
                          data = sample_B1,
                          verbose = TRUE,
                          svydesign = sample_A_svy_cal,
                          control_selection = controlSel(h=1, est_method_sel = "gee", nfolds = 5),
                          control_inference = controlInf(vars_selection = TRUE))
  
  est_ipw_mle <- nonprob(selection =  X_formula,
                         target = ~ Y_11 + Y_12 + Y_21 + Y_22,
                         data = sample_B1,
                         svydesign = sample_A_svy_cal,
                         verbose = TRUE,
                         control_inference = controlInf(vars_selection = TRUE),
                         control_selection = controlSel(nfolds = 5))
  # ipw B2 ---------------------------------------------------------------------
  est_ipw_h1_0_2 <- nonprob(selection = X_formula,
                            target = ~ Y_11,
                            data = sample_B2,
                            svydesign = sample_A_svy_cal,
                            verbose = TRUE,
                            control_selection = controlSel(h=1, est_method_sel = "gee", nfolds = 5),
                            control_inference = controlInf(vars_selection = TRUE))
  
  est_ipw_h1_1_2 <- nonprob(selection = X_formula,
                            target = ~ Y_12,
                            data = sample_B2,
                            svydesign = sample_A_svy_cal,
                            verbose = TRUE,
                            control_selection = controlSel(h=1, est_method_sel = "gee", nfolds = 5),
                            control_inference = controlInf(vars_selection = TRUE))
  
  est_ipw_h1_2_2 <- nonprob(selection = X_formula,
                            target = ~ Y_21,
                            data = sample_B2,
                            svydesign = sample_A_svy_cal,
                            verbose = TRUE,
                            control_selection = controlSel(h=1, est_method_sel = "gee", nfolds = 5),
                            control_inference = controlInf(vars_selection = TRUE))
  
  est_ipw_h1_3_2 <- nonprob(selection = X_formula,
                            target = ~ Y_22,
                            data = sample_B2,
                            svydesign = sample_A_svy_cal,
                            verbose = TRUE,
                            control_selection = controlSel(h=1, est_method_sel = "gee", nfolds = 5),
                            control_inference = controlInf(vars_selection = TRUE))
  
  est_ipw_mle_2 <- nonprob(selection = X_formula,
                           target = ~ Y_11 + Y_12 + Y_21 + Y_22,
                           data = sample_B2,
                           svydesign = sample_A_svy_cal,
                           verbose = TRUE,
                           control_inference = controlInf(vars_selection = TRUE),
                           control_selection = controlSel(nfolds = 5))
  # dr ----------------------------------------------------------------------
  dr_y11 <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B1) {
    if (x %in% c("Y_11", "Y_12")) {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                    selection = X_formula,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    family_outcome = "gaussian",
                    control_selection = controlSel(nfolds = 5),
                    control_outcome = controlOut(nfolds = 5),
                    verbose = TRUE,
                    control_inference = controlInf(vars_selection = TRUE))
    } else {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                    selection = X_formula,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    family_outcome = "binomial",
                    control_selection = controlSel(nfolds = 5),
                    control_outcome = controlOut(nfolds = 5),
                    verbose = TRUE,
                    control_inference = controlInf(vars_selection = TRUE))
    }
    data.frame(y=x, dr$output, dr$confidence_interval)
  })
  dr_y11_h1 <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B1) {
    if (x %in% c("Y_11", "Y_12")) {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                    selection = X_formula,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_selection = controlSel(est_method_sel = "gee", h = 1, nfolds = 5),
                    family_outcome = "gaussian",
                    control_outcome = controlOut(nfolds = 5),
                    verbose = TRUE,
                    control_inference = controlInf(vars_selection = TRUE))
    } else {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                    selection = X_formula,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_selection = controlSel(est_method_sel = "gee", h = 1, nfolds = 5),
                    family_outcome = "binomial",
                    control_outcome = controlOut(nfolds = 5),
                    verbose = TRUE,
                    control_inference = controlInf(vars_selection = TRUE))
    }
    data.frame(y=x, dr$output, dr$confidence_interval)
  })
  dr_y11_2 <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B2) {
    if (x %in% c("Y_11", "Y_12")) {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                    selection = X_formula,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_outcome = controlOut(nfolds = 5),
                    verbose = TRUE,
                    family_outcome = "gaussian",
                    control_inference = controlInf(vars_selection = TRUE))
    } else {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                    selection = X_formula,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_outcome = controlOut(nfolds = 5),
                    verbose = TRUE,
                    family_outcome = "binomial",
                    control_inference = controlInf(vars_selection = TRUE))
    }
    data.frame(y=x, dr$output, dr$confidence_interval)
  })
  dr_y11_h1_2 <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B2) {
    if (x %in% c("Y_11", "Y_12")) {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                    selection = X_formula,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_selection = controlSel(est_method_sel = "gee", h= 1, nfolds = 5),
                    family_outcome = "gaussian",
                    control_outcome = controlOut(nfolds = 5),
                    verbose = TRUE,
                    control_inference = controlInf(vars_selection = TRUE))
    } else {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ X3 + X4 + X5 + X6")),
                    selection = X_formula,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_selection = controlSel(est_method_sel = "gee", h=1, nfolds = 5),
                    family_outcome = "binomial",
                    control_outcome = controlOut(nfolds = 5),
                    verbose = TRUE,
                    control_inference = controlInf(vars_selection = TRUE))
    }
    data.frame(y=x, dr$output, dr$confidence_interval)
  })
  
  # dr (MM) ----------------------------------------------------------------------
  dr_y11_mm <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B1) {
    if (x %in% c("Y_11", "Y_12")) {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                    selection = X_formula,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    family_outcome = "gaussian",
                    control_outcome = controlOut(nfolds = 5),
                    verbose = TRUE,
                    control_inference = controlInf(bias_correction = TRUE, vars_selection = TRUE))
    } else {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                    selection = X_formula,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    family_outcome = "binomial",
                    control_outcome = controlOut(nfolds = 5),
                    verbose = TRUE,
                    control_inference = controlInf(bias_correction = TRUE, vars_selection = TRUE))
    }
    data.frame(y=x, dr$output, dr$confidence_interval)
  })
  dr_y11_h1_mm <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B1) {
    if (x %in% c("Y_11", "Y_12")) {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                    selection = X_formula,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_selection = controlSel(h = 1, nfolds = 5),
                    family_outcome = "gaussian",
                    control_outcome = controlOut(nfolds = 5),
                    verbose = TRUE,
                    control_inference = controlInf(bias_correction = TRUE, vars_selection = TRUE))
    } else {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                    selection = X_formula,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_selection = controlSel(h=1, nfolds = 5),
                    family_outcome = "binomial",
                    control_outcome = controlOut(nfolds = 5),
                    verbose = TRUE,
                    control_inference = controlInf(bias_correction = TRUE, vars_selection = TRUE))
    }
    data.frame(y=x, dr$output, dr$confidence_interval)
  })
  dr_y11_mm_2 <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B2) {
    if (x %in% c("Y_11", "Y_12")) {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                    selection = X_formula,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_selection = controlSel(h=2, nfolds = 5),
                    family_outcome = "gaussian",
                    control_outcome = controlOut(nfolds = 5),
                    verbose = TRUE,
                    control_inference = controlInf(bias_correction = TRUE, vars_selection = TRUE))
    } else {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                    selection = X_formula,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_selection = controlSel(h=2, nfolds = 5),
                    family_outcome = "binomial",
                    control_outcome = controlOut(nfolds = 5),
                    verbose = TRUE,
                    control_inference = controlInf(bias_correction = TRUE, vars_selection = TRUE))
    }
    data.frame(y=x, dr$output, dr$confidence_interval)
  })
  dr_y11_h1_mm_2 <- lapply(c("Y_11", "Y_12", "Y_21", "Y_22"), FUN = function(x, data=sample_B2) {
    if (x %in% c("Y_11", "Y_12")) {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                    selection = X_formula,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_selection = controlSel(h=1, nfolds = 5),
                    family_outcome = "gaussian",
                    control_outcome = controlOut(nfolds = 5),
                    verbose = TRUE,
                    control_inference = controlInf(bias_correction = TRUE, vars_selection = TRUE))
    } else {
      dr <- nonprob(outcome = as.formula(glue("{x} ~ {as.character(X_formula)[2]}")),
                    selection = X_formula,
                    data = data,
                    svydesign = sample_A_svy_cal,
                    control_selection = controlSel(h=1, nfolds = 5),
                    family_outcome = "binomial",
                    control_outcome = controlOut(nfolds = 5),
                    verbose = TRUE,
                    control_inference = controlInf(bias_correction = TRUE, vars_selection = TRUE))
    }
    data.frame(y=x, dr$output, dr$confidence_interval)
  })
  # save results ------------------------------------------------------------
  results_var_sel[[b]] <- rbind(
    rbindlist(mass_imp_glm)[, ":="(dataset="PSM I", est="glm")],
    rbindlist(mass_imp_glm2)[, ":="(dataset="PSM II", est="glm")],
    rbindlist(mass_imp_glm_nn)[, ":="(dataset="PSM I", est="nn")],
    rbindlist(mass_imp_glm_nn2)[, ":="(dataset="PSM II", est="nn")],
    ## PSM I
    data.frame(y=rownames(est_ipw_h1_0$output),
               est_ipw_h1_0$output,est_ipw_h1_0$confidence_interval, dataset= "PSM I", est = "ipw (h=1)"),
    data.frame(y=rownames(est_ipw_h1_1$output),
               est_ipw_h1_1$output,est_ipw_h1_1$confidence_interval, dataset= "PSM I", est = "ipw (h=1)"),
    data.frame(y=rownames(est_ipw_h1_2$output),
               est_ipw_h1_2$output,est_ipw_h1_2$confidence_interval, dataset= "PSM I", est = "ipw (h=1)"),
    data.frame(y=rownames(est_ipw_h1_3$output),
               est_ipw_h1_3$output,est_ipw_h1_3$confidence_interval, dataset= "PSM I", est = "ipw (h=1)"),
    data.frame(y=rownames(est_ipw_mle$output),
               est_ipw_mle$output,est_ipw_mle$confidence_interval, dataset= "PSM I", est = "ipw (h=2)"),
    ## PSM II
    data.frame(y=rownames(est_ipw_h1_0_2$output),
               est_ipw_h1_0_2$output,est_ipw_h1_0_2$confidence_interval, dataset= "PSM II", est = "ipw (h=1)"),
    data.frame(y=rownames(est_ipw_h1_1_2$output),
               est_ipw_h1_1_2$output,est_ipw_h1_1_2$confidence_interval, dataset= "PSM II", est = "ipw (h=1)"),
    data.frame(y=rownames(est_ipw_h1_2_2$output),
               est_ipw_h1_2_2$output,est_ipw_h1_2_2$confidence_interval, dataset= "PSM II", est = "ipw (h=1)"),
    data.frame(y=rownames(est_ipw_h1_3_2$output),
               est_ipw_h1_3_2$output,est_ipw_h1_3_2$confidence_interval, dataset= "PSM II", est = "ipw (h=1)"),
    data.frame(y=rownames(est_ipw_mle_2$output),
               est_ipw_mle_2$output,est_ipw_mle_2$confidence_interval, dataset= "PSM II", est = "ipw (h=2)"),
    ## dr
    rbindlist(dr_y11)[, ":="(dataset="PSM I", est="dr (h=2)")],
    rbindlist(dr_y11_h1)[, ":="(dataset="PSM I", est="dr (h=1)")],
    rbindlist(dr_y11_2)[, ":="(dataset="PSM II", est="dr (h=2)")],
    rbindlist(dr_y11_h1_2)[, ":="(dataset="PSM II", est="dr (h=1)")],
    ## dr (mm)
    rbindlist(dr_y11_mm)[, ":="(dataset="PSM I", est="dr mm (h=2)")],
    rbindlist(dr_y11_h1_mm)[, ":="(dataset="PSM I", est="dr mm (h=1)")],
    rbindlist(dr_y11_mm_2)[, ":="(dataset="PSM II", est="dr mm (h=2)")],
    rbindlist(dr_y11_h1_mm_2)[, ":="(dataset="PSM II", est="dr mm (h=1)")]
  )
}

## collect errors
results_var_sel_df <- rbindlist(results_var_sel, idcol = "b")
results_var_sel_df[y == "Y_11", true := mean(Y_11)]
results_var_sel_df[y == "Y_12", true := mean(Y_12)]
results_var_sel_df[y == "Y_21", true := mean(Y_21)]
results_var_sel_df[y == "Y_22", true := mean(Y_22)]
saveRDS(results_var_sel_df, file = "results/yang2020-sel-x-500.rds")

## coverage
results_var_sel_df[!is.na(lower_bound), 
                     .(m = mean(lower_bound < true & upper_bound > true)), keyby=.(y, dataset, est)]

ggplot(data = results_var_sel_df[!is.na(lower_bound) & lower_bound > 0], 
       aes(x = est, y = mean)) + 
  geom_boxplot(position = "dodge") +
  geom_hline(aes(yintercept = true), color = "red", linetype = "dashed") +
  facet_wrap(dataset~y, ncol = 4, scales = "free_y") +
  theme(axis.text.x=element_text(angle=45, vjust=1,hjust=1))


results_correct_x_df[y == "Y_11"
                     , .(bias = (mean(mean) - mean(true))/mean(true)*100,
                         var = var(mean),
                         rmse = sqrt( (mean(mean) - mean(true))^2 + var(mean))), 
                     keyby=.(y, dataset, est)]

