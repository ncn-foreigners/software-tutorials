library(nonprobsvy)
library(survey)
library(sampling)
library(glue)
library(data.table)
library(ggplot2)

# data --------------------------------------------------------------------

set.seed(123123123)
N <- 10000
n_a <- 500
n_b <- 1000
n_b1 <- 0.7*n_b
n_b2 <- 0.3*n_b
x1 <- rnorm(N, 2, 1)
x2 <- rnorm(N, 2, 1)
y1 <- rnorm(N, 0.3 + 2*x1+ 2*x2, 1)
y2 <- rnorm(N, 0.3 + 0.5*x1^2+ 0.5*x2^2, 1)
strata <- x1 <= 2
pop <- data.frame(x1, x2, y1, y2, strata)

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
  
  res_y_prob <- colMeans(sample_a[, c("y1", "y2")])
  res_y_non <- colMeans(sample_b[, c("y1", "y2")])
  
  res_y_glm <- nonprob(outcome = y1 + y2 ~ x1 + x2, data = sample_b, svydesign = svy_a) 
  res_y_npar <- nonprob(outcome = y1 + y2 ~ x1 + x2, 
                        data = sample_b, 
                        svydesign = svy_a, 
                        method_outcome = "npar")
  
  results[[b]] <- data.frame(y = rep(1:2, times = 4),
                             mean = c(res_y_glm$output$mean, res_y_npar$output$mean, 
                                      res_y_prob, res_y_non),
                             SE = c(res_y_glm$output$SE,res_y_npar$output$SE,
                                    rep(NA, 4)),
                             lower_bound = c(res_y_glm$confidence_interval$lower_bound,
                                             res_y_npar$confidence_interval$lower_bound,
                                             rep(NA, 4)),
                             upper_bound = c(res_y_glm$confidence_interval$upper_bound, 
                                             res_y_npar$confidence_interval$upper_bound,
                                             rep(NA, 4)),
                             est = rep(c("glm", "npar", "prob", "naive"), each=2))
  
}



## collect errors
results_paper <- rbindlist(results, idcol = "b")
results_paper[y == 1, true := mean(y1)]
results_paper[y == 2, true := mean(y2)]
saveRDS(results_paper, file = "results/chen2022-replicates.rds")
