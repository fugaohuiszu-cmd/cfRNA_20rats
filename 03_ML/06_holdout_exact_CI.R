# Exact 95% confidence intervals from the saved hold-out results.
options(stringsAsFactors = FALSE)

table_dir <- "F:/DE-analysis/final/cfRNA_project/03_results/03_ML/02_strict_model/03_model_evaluation/tables"
input_file <- file.path(table_dir, "Model_performance.csv")
output_file <- file.path(table_dir, "Holdout_exact_CI.csv")
confidence_level <- 0.95

performance <- read.csv(input_file, check.names = FALSE)
holdout <- performance[performance$Set == "Validation set", ]
stopifnot(
  nrow(holdout) == 10L, !anyDuplicated(holdout$Model),
  all(holdout$N == 8L),
  all(holdout$TP + holdout$FN == 4L),
  all(holdout$TN + holdout$FP == 4L)
)

results <- do.call(rbind, lapply(seq_len(nrow(holdout)), function(i) {
  one <- holdout[i, ]
  successes <- c(one$TP, one$TN, one$TP + one$TN)
  totals <- c(one$TP + one$FN, one$TN + one$FP, one$N)
  metrics <- c("Sensitivity", "Specificity", "Accuracy")
  intervals <- t(vapply(seq_along(metrics), function(j) {
    as.numeric(stats::binom.test(
      successes[j], totals[j], conf.level = confidence_level
    )$conf.int)
  }, numeric(2)))
  estimates <- successes / totals
  stopifnot(isTRUE(all.equal(
    estimates, as.numeric(unlist(one[1, metrics], use.names = FALSE)),
    tolerance = 1e-12
  )))
  data.frame(
    Set = "Internal hold-out set", Model = one$Model, Metric = metrics,
    Successes = successes, N = totals, Estimate = estimates,
    CI_low = intervals[, 1], CI_high = intervals[, 2]
  )
}))

write.csv(results, output_file, row.names = FALSE)
cat("Models: 10; hold-out samples: 8; intervals: ", nrow(results), "\n", sep = "")
print(results[results$Model == "PAM", ], row.names = FALSE)
cat("Output: ", output_file, "\n", sep = "")
