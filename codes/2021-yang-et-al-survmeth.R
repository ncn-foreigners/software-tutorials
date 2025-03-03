library(nonprobsvy)
library(survey)
library(sampling)
library(glue)
library(data.table)
library(ggplot2)

# data --------------------------------------------------------------------

set.seed(123123123)
N <- 1000000
n_a <- 1000
x1 <- rnorm(N,1,1)
x2 <- rexp(N,1)
alp <- rnorm(N)
epsilon <- rnorm(N)
y11 <- 1 + x1 + x2 + alp + epsilon
y12 <- 0.5*(x1-1.5)^2 + x2^2 + alp + epsilon
y21 <- rbinom(N,1,plogis(1 + x1 + x2 + alp))
y22 <- rbinom(N,1,plogis(0.5*(x1-1.5)^2 + x2^2 + alp))
p1 <- plogis(x2)
p2 <- plogis(-3+(x1-1.5)^2+(x2-2)^2)
pop <- data.frame(x1,x2,y11,y12,y21,y22,p1,p2) |> setDT()

y_vars <- c("y11", "y12", "y21", "y22")

B <- 1000
results <- list()

for (b in 1:B) {
  
  set.seed(b)
  print(glue("Iteration: {b} out of {B}."))
  sample_a <- pop[sample(1:N, n_a),]
  sample_a$w_a <- N/n_a
  sample_a_svy <- svydesign(ids=~1, weights=~w_a, data=sample_a)
  
  sample_b1 <- pop[rbinom(n = N, size = 1, prob = pop$p1) == 1, ]
  sample_b2 <- pop[rbinom(n = N, size = 1, prob = pop$p2) == 1, ]
  
  ## 
  ys_hat <- svymean(~y11+y12+y21+y22, sample_a_svy)
  ys_naive_b1 <- colMeans(sample_b1[, c("y11", "y12", "y21", "y22")])
  ys_naive_b2 <- colMeans(sample_b2[, c("y11", "y12", "y21", "y22")])
  
  ## NN estimators only
  ys_hat_nni_b1 <- nonprob(data=sample_b1, 
                           outcome=y11+y12+y21+y22~x1 + x2,
                           svydesign = sample_a_svy,
                           method_outcome = "nn",
                           control_outcome = control_out(k=1))
  ys_hat_kni_b1 <- nonprob(data=sample_b1, 
                           outcome=y11+y12+y21+y22~x1 + x2,
                           svydesign = sample_a_svy,
                           method_outcome = "nn",
                           control_outcome = control_out(k=5))
  ys_hat_nni_b2 <- nonprob(data=sample_b2, 
                           outcome=y11+y12+y21+y22~x1 + x2,
                           svydesign = sample_a_svy,
                           method_outcome = "nn",
                           control_outcome = control_out(k=1))
  ys_hat_kni_b2 <- nonprob(data=sample_b2, 
                           outcome=y11+y12+y21+y22~x1 + x2,
                           svydesign = sample_a_svy,
                           method_outcome = "nn",
                           control_outcome = control_out(k=5))
  
  results[[b]] <- data.frame(y = rep(y_vars, times = 2*8),
                             source = rep(c("b1", "b2"), each = 4, times = 4),
                             mean = c(as.numeric(ys_hat), 
                                      as.numeric(ys_naive_b1), 
                                      ys_hat_nni_b1$output$mean, 
                                      ys_hat_kni_b1$output$mean, 
                                      as.numeric(ys_hat),
                                      as.numeric(ys_naive_b2),
                                      ys_hat_nni_b2$output$mean,
                                      ys_hat_kni_b2$output$mean),
                          SE = c(SE(ys_hat), 
                                 rep(NA, 4), 
                                 ys_hat_nni_b1$output$SE, 
                                 ys_hat_kni_b1$output$SE, 
                                 SE(ys_hat),
                                 rep(NA, 4),
                                 ys_hat_nni_b2$output$SE,
                                 ys_hat_kni_b2$output$SE),
                          lower_bound = c(confint(ys_hat)[, 1], 
                                          rep(NA, 4), 
                                          ys_hat_nni_b1$confidence_interval$lower_bound, 
                                          ys_hat_kni_b1$confidence_interval$lower_bound, 
                                          confint(ys_hat)[, 1],
                                          rep(NA, 4),
                                          ys_hat_nni_b2$confidence_interval$lower_bound, 
                                          ys_hat_kni_b2$confidence_interval$lower_bound),
                          upper_bound = c(confint(ys_hat)[, 2], 
                                          rep(NA, 4), 
                                          ys_hat_nni_b1$confidence_interval$upper_bound, 
                                          ys_hat_kni_b1$confidence_interval$upper_bound, 
                                          confint(ys_hat)[, 2],
                                          rep(NA, 4),
                                          ys_hat_nni_b2$confidence_interval$upper_bound, 
                                          ys_hat_kni_b2$confidence_interval$upper_bound),
                          est = rep(c("ht", "naive", "nn", "knn"), each = 4, times = 4))
}

## collect errors
results_paper <- rbindlist(results, idcol = "b")
results_paper[y == "y11", true := mean(y11)]
results_paper[y == "y12", true := mean(y12)]
results_paper[y == "y21", true := mean(y21)]
results_paper[y == "y22", true := mean(y22)]

saveRDS(results_paper, file = "results/yang2021-replicates.rds")
