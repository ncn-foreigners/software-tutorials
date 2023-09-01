library(nonprobsvy)
library(survey)
library(sampling)
library(glue)
library(data.table)
library(ggplot2)

# data --------------------------------------------------------------------

set.seed(seed_for_sim)
N <- 100000
n_a <- 500
n_b <- 1000
n_b1 <- 0.7*n_b
n_b2 <- 0.3*n_b
x <- rnorm(N, 2, 1)
e <- rnorm(N)
y1 <- 1 + 2*x + e
y2 <- 3 + x + 2*e
y3 <- 2.5 + 0.5*x^2 + e
strata <- x <= 2
pop <- data.frame(x, y1, y2, y3, strata)

B <- 1000
results <- list()

for (b in 1:B) {
  set.seed(b)
  if (b %% 50 == 0) print(glue("Iteration: {b} out of {B}."))
  sample_a <- pop[sample(1:N, n_a),]
  sample_a$w_a <- N/n_a
  svy_a <- svydesign(ids= ~1, weights = ~ w_a, data = sample_a)
  pop1 <- subset(pop, strata == TRUE)
  pop2 <- subset(pop, strata == FALSE)
  sample_b <- rbind(pop1[sample(1:nrow(pop1), n_b1), ],
                    pop2[sample(1:nrow(pop2), n_b2), ])
  
  results_y <- lapply(c("y1", "y2", "y3"), function(X) {
    
    form_y_mo <- as.formula(glue("{X} ~ x"))
    form_y_ipw <- as.formula(glue("~ {X}"))
    
    res_y1_glm <- nonprob(outcome = form_y_mo, data = sample_b, svydesign = svy_a) ## from the paper
    res_y1_nn <- nonprob(outcome = form_y_mo, data = sample_b, svydesign = svy_a, 
                         method_outcome = "nn")
    res_y1_ipw1 <- nonprob(target = form_y_ipw,  selection = ~ x,  data = sample_b, svydesign = svy_a,
                           control_selection = controlSel(est_method_sel = "mle")) ## from the paper
    res_y1_ipw2 <- nonprob(target = form_y_ipw,  selection = ~ x,  data = sample_b, svydesign = svy_a,
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
    res_y1_dr1 <- nonprob(outcome = form_y_mo,  selection = ~ x,  data = sample_b, svydesign = svy_a,
                          control_selection = controlSel(est_method_sel = "mle"))
    res_y1_dr2 <- nonprob(outcome = form_y_mo,  selection = ~ x,  data = sample_b, svydesign = svy_a,
                          control_selection = controlSel(est_method_sel = "gee", h = 1))
    res_y1_dr3 <- nonprob(outcome = form_y_mo,  selection = ~ x,  data = sample_b, svydesign = svy_a,
                          control_selection = controlSel(est_method_sel = "mle"),
                          control_inference = controlInf(bias_correction = T))
    res_y1_dr4 <- nonprob(outcome = form_y_mo,  selection = ~ x,  data = sample_b, svydesign = svy_a,
                          control_selection = controlSel(est_method_sel = "gee", h = 1),
                          control_inference = controlInf(bias_correction = T))
    
    prob_res <- svymean(~y1, svy_a)
    prob_res_df <- as.data.frame(prob_res)
    prob_res_ci <- confint(prob_res)
    
    rbind(data.frame(mean=mean(sample_b[,X]), SE=NA, lower_bound=NA, upper_bound=NA,est='naive'),
          data.frame(mean=prob_res_df[,1],SE=prob_res_df[,2], lower_bound=prob_res_ci[1], upper_bound=prob_res_ci[2], est='prob'),
          cbind(res_y1_glm$output, res_y1_glm$confidence_interval, est='mi (glm)'),
          cbind(res_y1_nn$output, res_y1_nn$confidence_interval, est='mi (nn)'),
          cbind(res_y1_ipw1$output, res_y1_ipw1$confidence_interval, est='ipw (mle)'),
          cbind(res_y1_ipw2$output, res_y1_ipw2$confidence_interval, est='gee (h=1)'),
          cbind(res_y1_dr1$output, res_y1_dr1$confidence_interval, est='dr (mle)'),
          cbind(res_y1_dr2$output, res_y1_dr2$confidence_interval, est='dr (gee,h=1)'),
          cbind(res_y1_dr3$output, res_y1_dr3$confidence_interval, est='dr (mle,bias)'),
          cbind(res_y1_dr4$output, res_y1_dr4$confidence_interval, est='dr (gee,h=1,bias)'))
  })
  
  results[[b]] <- rbindlist(results_y, idcol = "y")
  
}



## collect errors
results_paper <- rbindlist(results, idcol = "b")
results_paper[y == 1, true := mean(y1)]
results_paper[y == 2, true := mean(y2)]
results_paper[y == 3, true := mean(y3)]
saveRDS(results_paper, file = "results/kim2021-replicates.rds")

## coverage
results_paper[!is.na(lower_bound), .(m = mean(lower_bound < true & upper_bound > true)), keyby=.(y, est)]

## plot
ggplot(data = results_paper, aes(x = est, y = mean)) + 
  geom_boxplot(position = "dodge") +
  geom_hline(aes(yintercept = true), color = "red", linetype = "dashed") +
  facet_wrap(~y, ncol = 3, scales = "free_y") +
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1))


## bias, var and mse
results_paper[, .(bias = mean(mean) - mean(true), var = var(mean)*1000,
                  mse = (mean(mean) - mean(true))^2 + var(mean),
                  rmse = sqrt((mean(mean) - mean(true))^2 + var(mean))), 
              keyby=.(y, est)][, ReMSE:=mse/mse[est=="prob"]*100][]
