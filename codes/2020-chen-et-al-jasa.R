library(nonprobsvy)
library(survey)
library(sampling)
library(glue)
library(data.table)
library(ggplot2)
library(Rcpp)

## based on https://sas.uwaterloo.ca/~cbwu/Rcodes/SystematicPPS.txt (returns 0/1 instead of numbers)

cppFunction('
IntegerVector syspps_rcpp_indicator(NumericVector x, int n) {
  int N = x.size();
  double total_x = sum(x);
  
  // Create indices 1 to N
  IntegerVector U(N);
  for (int i = 0; i < N; i++) {
    U[i] = i + 1;
  }
  
  // Shuffle indices
  U = sample(U, N, false, R_NilValue);
  
  // Reorder x values (adjusting for 1-based indices)
  NumericVector xx(N);
  for (int i = 0; i < N; i++) {
    xx[i] = x[U[i] - 1];
  }
  
  // Calculate cumulative sum
  NumericVector cum_xx(N);
  cum_xx[0] = xx[0];
  for (int i = 1; i < N; i++) {
    cum_xx[i] = cum_xx[i-1] + xx[i];
  }
  
  // Calculate z values
  NumericVector z(N);
  for (int i = 0; i < N; i++) {
    z[i] = n * cum_xx[i] / total_x;
  }
  
  // Select sample
  double r = R::runif(0, 1);
  IntegerVector selected_indices(n);
  int count = 0;
  
  for (int i = 0; i < N && count < n; i++) {
    if (z[i] >= r) {
      selected_indices[count] = U[i];
      count++;
      r += 1.0;
    }
  }
  
  // Create indicator vector (0/1)
  IntegerVector indicator(N, 0);  // Initialize all elements to 0
  
  // Set selected elements to 1
  for (int i = 0; i < count; i++) {
    // Adjust for 0-based indexing in C++
    indicator[selected_indices[i] - 1] = 1;
  }
  
  return indicator;
}')

set.seed(123)
# sizes of population and probability sample
N <- 20000 # population
n_a <- 500 # nonprob
n_b <- 1000 # probability
# data
z1 <- rbinom(N, 1, 0.7)
z2 <- runif(N, 0, 2)
z3 <- rexp(N, 1)
z4 <- rchisq(N, 4)

# covariates
x1 <- z1
x2 <- z2 + 0.3 * z2
x3 <- z3 + 0.2 * (z1 + z2)
x4 <- z4 + 0.1 * (z1 + z2 + z3)
epsilon <- rnorm(N)
sigma_30 <- 10.4
sigma_50 <- 5.2
sigma_80 <- 2.4

# response variables
y30 <- 2 + x1 + x2 + x3 + x4 + sigma_30 * epsilon
y50 <- 2 + x1 + x2 + x3 + x4 + sigma_50 * epsilon
y80 <- 2 + x1 + x2 + x3 + x4 + sigma_80 * epsilon

# population
bb <- uniroot(f = function(x) sum(plogis(x + 0.1 * x1 + 0.2 * x2 + 0.1 * x3 + 0.2 * x4)) - 500, 
              lower=-10, upper=10)

sim_data <- data.frame(y30, y50, y80, x1, x2, x3, x4,
                       rho = plogis(bb$root + 0.1 * x1 + 0.2 * x2 + 0.1 * x3 + 0.2 * x4),
                       p_prob=inclusionprobabilities(x3 + 0.2051, n = n_a))

ys_names <- paste0("y",c("30","50","80"))

ipw_controls <- control_sel(est_method = "gee", 
                            nleqslv_method = "Newton",
                            nleqslv_global = "cline")


# sampling
R_sims <- 5000
res_sims <- list()

for (r in 1:R_sims) {
  # Progress tracking
  if (r %% 100 == 0) print(glue("Iteration: {r} out of {R_sims}."))
  
  # Try the entire iteration with error handling
  tryCatch({
    set.seed(r)
    
    # Sampling
    sim_data$flag_nonprob <- UPpoisson(sim_data$rho) ## sampling nonprob
    sim_data$flag_prob <- syspps_rcpp_indicator(sim_data$p_prob, n_b) ## sampling prob
    
    # Prepare samples
    sample_nonp <- subset(sim_data, flag_nonprob == 1)
    sample_prob <- subset(sim_data, flag_prob == 1)
    
    # Create survey design
    sample_prob_svy <- svydesign(ids=~1, probs=~p_prob, data=sample_prob)
    
    # Define models
    form1_ipw <- ~ x1 + x2 + x3 + x4
    form1_mi <- y30 + y50 + y80 ~ x1 + x2 + x3 + x4
    form2_ipw <- ~ x1 + x2 + x3
    form2_mi <- y30 + y50 + y80 ~ x1 + x2 + x3
    y_target <- ~y30 + y50 + y80
    
    # Initialize storage for this iteration's results
    results_list <- list()
    
    # Naive estimator - simple mean from nonprob sample
    results_list$naive <- tryCatch({
      colMeans(sample_nonp[, c("y30", "y50", "y80")])
    }, error = function(e) {
      print(glue("Iteration {r}: Error in naive estimator - {conditionMessage(e)}"))
      return(NULL)
    })
    
    # Probability sample estimator
    prob_results <- tryCatch({
      est_prob <- svymean(y_target, sample_prob_svy)
      est_prob_ci <- confint(est_prob)
      list(est = est_prob, ci = est_prob_ci)
    }, error = function(e) {
      print(glue("Iteration {r}: Error in probability estimator - {conditionMessage(e)}"))
      return(NULL)
    })
    results_list$prob <- prob_results
    
    # MI with correct specification
    results_list$mi_cor <- tryCatch({
      nonprob(outcome = form1_mi, data = sample_nonp, svydesign = sample_prob_svy)
    }, error = function(e) {
      print(glue("Iteration {r}: Error in MI correct - {conditionMessage(e)}"))
      return(NULL)
    })
    
    # MI with incorrect specification
    results_list$mi_inc <- tryCatch({
      nonprob(outcome = form2_mi, data = sample_nonp, svydesign = sample_prob_svy)
    }, error = function(e) {
      print(glue("Iteration {r}: Error in MI incorrect - {conditionMessage(e)}"))
      return(NULL)
    })
    
    # IPW with correct specification
    results_list$ipw_cor <- tryCatch({
      nonprob(selection = form1_ipw, target = y_target, data = sample_nonp, 
              svydesign = sample_prob_svy, control_selection = ipw_controls)
    }, error = function(e) {
      print(glue("Iteration {r}: Error in IPW correct - {conditionMessage(e)}"))
      return(NULL)
    })
    
    # IPW with incorrect specification
    results_list$ipw_inc <- tryCatch({
      nonprob(selection = form2_ipw, target = y_target, data = sample_nonp, 
              svydesign = sample_prob_svy, control_selection = ipw_controls)
    }, error = function(e) {
      print(glue("Iteration {r}: Error in IPW incorrect - {conditionMessage(e)}"))
      return(NULL)
    })
    
    # DR with both models correctly specified (TT)
    results_list$dr_tt <- tryCatch({
      nonprob(outcome = form1_mi, selection = form1_ipw, target = y_target, 
              data = sample_nonp, svydesign = sample_prob_svy, 
              control_selection = ipw_controls)
    }, error = function(e) {
      print(glue("Iteration {r}: Error in DR (TT) - {conditionMessage(e)}"))
      return(NULL)
    })
    
    # DR with outcome incorrect, selection correct (FT)
    results_list$dr_ft <- tryCatch({
      nonprob(outcome = form2_mi, selection = form1_ipw, target = y_target, 
              data = sample_nonp, svydesign = sample_prob_svy, 
              control_selection = ipw_controls)
    }, error = function(e) {
      print(glue("Iteration {r}: Error in DR (FT) - {conditionMessage(e)}"))
      return(NULL)
    })
    
    # DR with outcome correct, selection incorrect (TF)
    results_list$dr_tf <- tryCatch({
      nonprob(outcome = form1_mi, selection = form2_ipw, target = y_target, 
              data = sample_nonp, svydesign = sample_prob_svy, 
              control_selection = ipw_controls)
    }, error = function(e) {
      print(glue("Iteration {r}: Error in DR (TF) - {conditionMessage(e)}"))
      return(NULL)
    })
    
    # DR with both incorrectly specified (FF)
    results_list$dr_ff <- tryCatch({
      nonprob(outcome = form2_mi, selection = form2_ipw, target = y_target, 
              data = sample_nonp, svydesign = sample_prob_svy, 
              control_selection = ipw_controls)
    }, error = function(e) {
      print(glue("Iteration {r}: Error in DR (FF) - {conditionMessage(e)}"))
      return(NULL)
    })
    
    # Build results dataframe only with successful models
    result_rows <- list()
    
    # Add naive estimator results if available
    if (!is.null(results_list$naive)) {
      result_rows$naive <- cbind(
        y = ys_names, 
        data.frame(
          mean = results_list$naive, 
          SE = NA, 
          lower_bound = NA, 
          upper_bound = NA
        ), 
        est = "naive"
      )
    }
    
    # Add probability estimator results if available
    if (!is.null(results_list$prob)) {
      result_rows$prob <- cbind(
        y = ys_names, 
        data.frame(
          mean = as.vector(results_list$prob$est), 
          SE = unname(SE(results_list$prob$est)), 
          lower_bound = results_list$prob$ci[,1], 
          upper_bound = results_list$prob$ci[,2]
        ), 
        est = "prob"
      )
    }
    
    # Function to safely add nonprob results to the result_rows list
    add_nonprob_results <- function(result, name) {
      if (!is.null(result)) {
        result_rows[[name]] <<- cbind(
          y = ys_names, 
          result$output, 
          result$confidence_interval, 
          est = name
        )
      }
    }
    
    # Add all nonprob model results if available
    add_nonprob_results(results_list$mi_cor, "mi (T)")
    add_nonprob_results(results_list$mi_inc, "mi (F)")
    add_nonprob_results(results_list$ipw_cor, "ipw (T)")
    add_nonprob_results(results_list$ipw_inc, "ipw (F)")
    add_nonprob_results(results_list$dr_tt, "dr (TT)")
    add_nonprob_results(results_list$dr_ft, "dr (FT)")
    add_nonprob_results(results_list$dr_tf, "dr (TF)")
    add_nonprob_results(results_list$dr_ff, "dr (FF)")
    
    # Combine all result rows if any exist
    if (length(result_rows) > 0) {
      res_sims[[r]] <- do.call(rbind, result_rows)
    } else {
      print(glue("Iteration {r}: No successful models"))
      res_sims[[r]] <- NULL
    }
    
  }, error = function(e) {
    # Catch any other errors in the iteration
    print(glue("Iteration {r}: Fatal error - {conditionMessage(e)}"))
    res_sims[[r]] <- NULL
  })
}

# Clean up the results list by removing NULL entries
res_sims <- res_sims[!sapply(res_sims, is.null)]
results_paper <- rbindlist(res_sims, idcol = "r")
results_paper <- results_paper[!r %in% unique(results_paper$r[is.na(results_paper$mean)])]

# Summary of successful iterations
print(glue("Successfully completed {length(unique(results_paper$r))} out of {R_sims} iterations"))

## collect errors

results_paper[y == "y30", true := mean(y30)]
results_paper[y == "y50", true := mean(y50)]
results_paper[y == "y80", true := mean(y80)]

saveRDS(results_paper, file = "results/chen2020-replicates.rds")


