
# ============================================================================
# Copyright (c) 2025 Thomas Broggini and Liu Xiao. All rights reserved.
#
# This code was jointly written by:
# Thomas Broggini (Frankfurt University, Germany)
# Liu Xiao (Xiangyang First Hospital, China)
#
# The way how the mathematical models and methodologies implemented in this code have been
# individually customized and are not intended for generic use. 
#
# For permission requests or inquiries, please contact:
# Thomas Broggini  : broggini@med.uni-frankfurt.de
# Liu Xiao         : careyneurosurgery@gmail.com  /  carey-lau@foxmail.com
#
# Unauthorized use will be considered a violation of intellectual property rights.
# ============================================================================



#!/usr/bin/env Rscript

library(poweRlaw)
library(jsonlite)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
input_csv <- args[1]
output_json <- args[2]

# read data from matlab
data <- read.csv(input_csv)
k <- data$k
freq <- data$freq
data_vector <- as.numeric(rep(k, freq))
if (any(is.na(data_vector))) {
  stop("Non-numeric values detected in degree data. Check input CSV.")
}

# pre-powerlaw fitting，find optimal xmin
m <- displ$new(data_vector)
m$setXmin(2)  # begain xmin=2
est <- estimate_xmin(m)
xmin_initial <- est$xmin
alpha_initial <- est$pars

# make sure xmin has >=50 data points，or find xmin meet 50 points 
n_points_above_xmin <- sum(data_vector >= xmin_initial)
unique_k <- sort(unique(data_vector))

xmin_adjusted <- xmin_initial
if (n_points_above_xmin < 50) {
 
  k_candidate <- xmin_initial
  found <- FALSE
  while (k_candidate >= min(unique_k) && !found) {
    if (sum(data_vector >= k_candidate) >= 50) {
      xmin_adjusted <- k_candidate
      found <- TRUE
    }
    k_candidate <- k_candidate - 1
  }
  
  if (!found) {
    k_candidate <- xmin_initial + 1
    while (k_candidate <= max(unique_k) && !found) {
      if (sum(data_vector >= k_candidate) >= 50) {
        xmin_adjusted <- k_candidate
        found <- TRUE
      }
      k_candidate <- k_candidate + 1
    }
  }
  # if cannot find xmin has 50 points, warning 
  if (!found) {
    xmin_adjusted <- xmin_initial
    warning("No valid xmin found with >=50 data points. Using initial xmin.")
  }
  # update
  m$setXmin(xmin_adjusted)
  est_pars <- estimate_pars(m)
  m$setPars(est_pars)
} else {
  # 
  m$setXmin(xmin_initial)
  m$setPars(alpha_initial)
}

alpha_adjusted <- m$pars
data_cut <- data_vector[data_vector >= xmin_adjusted]

# bootstrap KS，2000，15
p1_result <- bootstrap_p(m, no_of_sims = 2000, threads = 15)
p1 <- p1_result$p

# possion vs powerlaw
m_pois <- dispois$new(data_vector)
m_pois$setXmin(xmin_adjusted)
est_pois_pars <- estimate_pars(m_pois)
m_pois$setPars(est_pois_pars)

# compare_distributions
comp_pois <- compare_distributions(m, m_pois)

LR <- comp_pois$test_statistic
p_value_lr <- comp_pois$p_two_sided

# trunca powerlaw，>=20
if (sum(data_vector >= xmin_adjusted) >= 20) {
  m_trunc <- disexp$new(data_vector)
  m_trunc$setXmin(xmin_adjusted)
  est_trunc_pars <- estimate_pars(m_trunc)
  m_trunc$setPars(est_trunc_pars)
  
  p1_trunc <- bootstrap_p(m_trunc, no_of_sims = 1000, threads = 15)$p
  
  comp_trunc <- compare_distributions(m, m_trunc)
  LR_trunc <- comp_trunc$test_statistic
  p_value_lr_trunc <- comp_trunc$p_two_sided
} else {
  m_trunc <- NULL
  p1_trunc <- NA
  p_value_lr_trunc <- NA
  LR_trunc <- NA
}

# output
result <- list(
  xmin_initial = xmin_initial,
  xmin_adjusted = xmin_adjusted,
  alpha_initial = alpha_initial,
  alpha_adjusted = alpha_adjusted,
  p_value_ks = p1,
  p_value_lr = p_value_lr,
  powerlaw_plausible = (p1 > 0.1),
  powerlaw_better = (p_value_lr < 0.025),
  truncated_alpha = if (!is.null(m_trunc)) m_trunc$pars[1] else NA,
  truncated_lambda = if (!is.null(m_trunc)) m_trunc$pars[2] else NA,
  p_value_ks_trunc = p1_trunc,
  p_value_lr_trunc = p_value_lr_trunc,
  trunc_better_than_power = if (!is.null(m_trunc)) (p_value_lr_trunc < 0.05) else NA
)

write_json(result, output_json, auto_unbox = TRUE)
