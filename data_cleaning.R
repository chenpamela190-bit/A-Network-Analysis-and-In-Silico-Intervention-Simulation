# Data Cleaning Script for SCL-90-R Mental Health Screening Data
# Implements two-stage data quality filtering protocol
# Author: Generated based on research methods
# Date: 2026-02-11

# Load required libraries
library(careless)  # For longstring analysis
library(psych)     # For Mahalanobis distance calculation

# ============================================================================
# STAGE 1: LONGSTRING ANALYSIS
# ============================================================================
# Purpose: Identify and exclude participants showing non-engagement through
# consecutive identical responses (indicating careless responding)
# Exclusion criterion: longstring >= 15
# ============================================================================

#' Calculate Longstring for Each Participant
#'
#' @param data A data frame containing item responses
#' @param item_cols Vector of column names or indices for SCL-90-R items
#' @return Vector of maximum consecutive identical responses per participant
calculate_longstring <- function(data, item_cols) {
  # Use careless package's longstring function
  longstring_values <- longstring(data[, item_cols])
  return(longstring_values)
}

#' Flag Participants Based on Longstring Threshold
#'
#' @param data A data frame containing item responses
#' @param item_cols Vector of column names or indices for SCL-90-R items
#' @param threshold Maximum acceptable longstring value (default: 15)
#' @return List containing flagged indices, clean data, and exclusion summary
stage1_longstring_filter <- function(data, item_cols, threshold = 15) {

  cat("=== STAGE 1: LONGSTRING ANALYSIS ===\n")
  cat("Calculating maximum consecutive identical responses...\n")

  # Calculate longstring values
  data$longstring <- calculate_longstring(data, item_cols)

  # Identify cases exceeding threshold
  flagged_indices <- which(data$longstring >= threshold)
  n_flagged <- length(flagged_indices)
  pct_flagged <- (n_flagged / nrow(data)) * 100

  # Create clean dataset
  data_clean <- data[data$longstring < threshold, ]

  # Summary statistics
  cat(sprintf("Initial sample size: N = %d\n", nrow(data)))
  cat(sprintf("Longstring threshold: >= %d\n", threshold))
  cat(sprintf("Cases flagged: n = %d (%.1f%%)\n", n_flagged, pct_flagged))
  cat(sprintf("Remaining sample: N = %d\n", nrow(data_clean)))
  cat(sprintf("Longstring range in clean data: %d to %d\n",
              min(data_clean$longstring), max(data_clean$longstring)))
  cat("\n")

  # Return results
  return(list(
    data_clean = data_clean,
    flagged_indices = flagged_indices,
    longstring_values = data$longstring,
    n_excluded = n_flagged,
    pct_excluded = pct_flagged,
    threshold = threshold
  ))
}


# ============================================================================
# STAGE 2: MULTIVARIATE OUTLIER DETECTION (MAHALANOBIS DISTANCE)
# ============================================================================
# Purpose: Identify extreme multivariate outliers whose response patterns
# deviate substantially from the sample centroid in 35-dimensional space
# Exclusion criterion: D² exceeding critical χ² value (df = 35, α = .001)
# ============================================================================

#' Calculate Mahalanobis Distance for Response Patterns
#'
#' @param data A data frame containing item responses
#' @param item_cols Vector of column names or indices for SCL-90-R items
#' @return Vector of Mahalanobis D² values for each participant
calculate_mahalanobis <- function(data, item_cols) {

  # Extract item data
  item_data <- data[, item_cols]

  # Calculate center (means) and covariance matrix
  center <- colMeans(item_data, na.rm = TRUE)
  cov_matrix <- cov(item_data, use = "pairwise.complete.obs")

  # Calculate Mahalanobis distance using psych package
  # Note: outlier() function from psych returns D² values
  mahal_d2 <- mahalanobis(item_data, center = center, cov = cov_matrix)

  return(mahal_d2)
}

#' Flag Multivariate Outliers Based on Mahalanobis Distance
#'
#' @param data A data frame containing item responses (output from Stage 1)
#' @param item_cols Vector of column names or indices for SCL-90-R items
#' @param df Degrees of freedom (number of items, default: 35)
#' @param alpha Significance level for chi-square cutoff (default: 0.001)
#' @return List containing flagged indices, clean data, and exclusion summary
stage2_mahalanobis_filter <- function(data, item_cols, df = 35, alpha = 0.001) {

  cat("=== STAGE 2: MULTIVARIATE OUTLIER DETECTION ===\n")
  cat("Calculating Mahalanobis distance (D²)...\n")

  # Calculate Mahalanobis D² values
  data$mahal_d2 <- calculate_mahalanobis(data, item_cols)

  # Determine critical value from chi-square distribution
  # D² follows χ² distribution under multivariate normality
  critical_value <- qchisq(p = 1 - alpha, df = df)

  # Identify cases exceeding critical value
  flagged_indices <- which(data$mahal_d2 > critical_value)
  n_flagged <- length(flagged_indices)
  pct_flagged <- (n_flagged / nrow(data)) * 100

  # Create clean dataset
  data_clean <- data[data$mahal_d2 <= critical_value, ]

  # Summary statistics
  cat(sprintf("Sample size entering Stage 2: N = %d\n", nrow(data)))
  cat(sprintf("Degrees of freedom: df = %d\n", df))
  cat(sprintf("Significance level: α = %.3f\n", alpha))
  cat(sprintf("Critical χ² value: %.2f\n", critical_value))
  cat(sprintf("Cases flagged: n = %d (%.1f%%)\n", n_flagged, pct_flagged))
  cat(sprintf("Final clean sample: N = %d\n", nrow(data_clean)))
  cat(sprintf("D² range in clean data: %.2f to %.2f\n",
              min(data_clean$mahal_d2), max(data_clean$mahal_d2)))
  cat(sprintf("D² mean (SD) in clean data: %.2f (%.2f)\n",
              mean(data_clean$mahal_d2), sd(data_clean$mahal_d2)))
  cat("\n")

  # Return results
  return(list(
    data_clean = data_clean,
    flagged_indices = flagged_indices,
    mahal_d2_values = data$mahal_d2,
    critical_value = critical_value,
    n_excluded = n_flagged,
    pct_excluded = pct_flagged,
    df = df,
    alpha = alpha
  ))
}


# ============================================================================
# COMBINED WORKFLOW: TWO-STAGE SEQUENTIAL FILTERING
# ============================================================================

#' Execute Complete Two-Stage Data Quality Filtering
#'
#' @param data A data frame containing raw SCL-90-R responses
#' @param item_cols Vector of column names or indices for the 35 selected items
#' @param longstring_threshold Threshold for Stage 1 (default: 15)
#' @param mahal_df Degrees of freedom for Stage 2 (default: 35)
#' @param mahal_alpha Significance level for Stage 2 (default: 0.001)
#' @return List containing clean data and comprehensive filtering summary
execute_data_cleaning <- function(data,
                                   item_cols,
                                   longstring_threshold = 15,
                                   mahal_df = 35,
                                   mahal_alpha = 0.001) {

  cat("╔════════════════════════════════════════════════════════════╗\n")
  cat("║  TWO-STAGE DATA QUALITY FILTERING PROTOCOL                ║\n")
  cat("║  Sequential filtering for careless responding detection   ║\n")
  cat("╚════════════════════════════════════════════════════════════╝\n\n")

  initial_n <- nrow(data)

  # STAGE 1: Longstring Analysis
  stage1_results <- stage1_longstring_filter(
    data = data,
    item_cols = item_cols,
    threshold = longstring_threshold
  )

  # STAGE 2: Mahalanobis Distance (applied to Stage 1 clean data)
  stage2_results <- stage2_mahalanobis_filter(
    data = stage1_results$data_clean,
    item_cols = item_cols,
    df = mahal_df,
    alpha = mahal_alpha
  )

  # Final summary
  final_n <- nrow(stage2_results$data_clean)
  total_excluded <- initial_n - final_n
  retention_rate <- (final_n / initial_n) * 100

  cat("=== FINAL SUMMARY ===\n")
  cat(sprintf("Initial sample: N = %d\n", initial_n))
  cat(sprintf("Stage 1 exclusions (longstring): n = %d (%.1f%%)\n",
              stage1_results$n_excluded, stage1_results$pct_excluded))
  cat(sprintf("Stage 2 exclusions (Mahalanobis): n = %d (%.1f%%)\n",
              stage2_results$n_excluded, stage2_results$pct_excluded))
  cat(sprintf("Total exclusions: n = %d (%.1f%%)\n",
              total_excluded, (total_excluded/initial_n)*100))
  cat(sprintf("Final analytic sample: N = %d (%.1f%% retention)\n",
              final_n, retention_rate))
  cat("══════════════════════════════════════════════════════════\n\n")

  # Compile comprehensive results
  results <- list(
    # Clean data
    data_clean = stage2_results$data_clean,

    # Stage-specific results
    stage1 = stage1_results,
    stage2 = stage2_results,

    # Overall summary
    summary = list(
      initial_n = initial_n,
      final_n = final_n,
      total_excluded = total_excluded,
      retention_rate = retention_rate,
      stage1_excluded = stage1_results$n_excluded,
      stage2_excluded = stage2_results$n_excluded
    )
  )

  return(results)
}


# ============================================================================
# EXAMPLE USAGE
# ============================================================================
#
# # Load your data
# raw_data <- read.csv("scl90r_data.csv")
#
# # Define the 35 selected items (adjust column names as needed)
# # Based on methods: 12 Somatization + 13 Depression + 10 Anxiety items
# selected_items <- c(
#   # Somatization (12 items)
#   "Q1", "Q4", "Q12", "Q27", "Q40", "Q42", "Q48", "Q49", "Q52", "Q53", "Q56", "Q58",
#   # Depression (13 items)
#   "Q5", "Q14", "Q15", "Q20", "Q22", "Q26", "Q29", "Q30", "Q31", "Q32", "Q54", "Q71", "Q79",
#   # Anxiety (10 items)
#   "Q2", "Q17", "Q23", "Q33", "Q39", "Q57", "Q72", "Q78", "Q80", "Q86"
# )
#
# # Execute two-stage filtering
# cleaning_results <- execute_data_cleaning(
#   data = raw_data,
#   item_cols = selected_items,
#   longstring_threshold = 15,
#   mahal_df = 35,
#   mahal_alpha = 0.001
# )
#
# # Extract clean data
# clean_data <- cleaning_results$data_clean
#
# # Save clean data
# write.csv(clean_data, "scl90r_clean_data.csv", row.names = FALSE)
#
# # Save filtering report
# saveRDS(cleaning_results, "data_cleaning_results.rds")
