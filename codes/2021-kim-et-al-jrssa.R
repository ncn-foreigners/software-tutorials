library(nonprobsvy)
library(survey)
library(sampling)
library(glue)
library(data.table)
library(ggplot2)

# data --------------------------------------------------------------------

set.seed(123123123)
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
  print(glue("Iteration: {b} out of {B}."))
  sample_a <- pop[sample(1:N, n_a),]
  sample_a$w_a <- N/n_a
  svy_a <- svydesign(ids= ~1, weights = ~ w_a, data = sample_a)
  pop1 <- subset(pop, strata == TRUE)
  pop2 <- subset(pop, strata == FALSE)
  sample_b <- rbind(pop1[sample(1:nrow(pop1), n_b1), ],
                    pop2[sample(1:nrow(pop2), n_b2), ])
  
  res_y_prob <- colMeans(sample_a[, c("y1", "y2", "y3")])
  res_y_non <- colMeans(sample_b[, c("y1", "y2", "y3")])
  
  res_y_glm <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_b, svydesign = svy_a) 
  
  res_y_ipw <- nonprob(target = ~y1 + y2 + y3,  
                       selection = ~ x,  data = sample_b, 
                       svydesign = svy_a,
                       control_selection = control_sel(est_method =  "mle")) 
  
  
  results[[b]] <- data.frame(y = rep(1:3, times = 4),
                             mean = c(res_y_glm$output$mean, res_y_ipw$output$mean, 
                                      res_y_prob, res_y_non),
                             SE = c(res_y_glm$output$SE,res_y_ipw$output$SE,
                                    rep(NA, 6)),
                             lower_bound = c(res_y_glm$confidence_interval$lower_bound, res_y_ipw$confidence_interval$lower_bound,
                                             rep(NA, 6)),
                             upper_bound = c(res_y_glm$confidence_interval$upper_bound, res_y_ipw$confidence_interval$upper_bound,
                                             rep(NA, 6)),
                             est = rep(c("mi", "ipw", "prob", "naive"), each=3))
  
}



## collect errors
results_paper <- rbindlist(results, idcol = "b")
results_paper[y == 1, true := mean(y1)]
results_paper[y == 2, true := mean(y2)]
results_paper[y == 3, true := mean(y3)]
saveRDS(results_paper, file = "results/kim2021-replicates.rds")