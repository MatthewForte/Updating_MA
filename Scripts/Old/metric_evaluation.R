# Monotone-Adjusted Jumpiness (MAJ)
monotone_adjusted_jumpiness <- function(x, 
                                        residual_jump_fun = c("MSSD","MASD","TotalVariation"),
                                        normalize_by = c("sd","mad","range","none")) {
  # x: numeric vector (length >= 2)
  # residual_jump_fun: how to measure jumpiness on residuals
  # normalize_by: how to normalize residual jumpiness to be scale-free
  # returns a list with M (isotonic R2), J (normalized residual jumpiness), MAJ = (1-M)*J
  
  if (length(x) < 2) stop("x must have length >= 2")
  residual_jump_fun <- match.arg(residual_jump_fun)
  normalize_by <- match.arg(normalize_by)
  
  n <- length(x)
  idx <- seq_len(n)
  
  # isotonic fits: non-decreasing and non-increasing (choose better)
  fit_inc <- stats::isoreg(idx, x)$yf
  ssres_inc <- sum((x - fit_inc)^2)
  sstot <- sum((x - mean(x))^2)
  r2_inc <- if (sstot > eps) 1 - ssres_inc/sstot else NA
  
  fit_dec <- -stats::isoreg(idx, -x)$yf
  ssres_dec <- sum((x - fit_dec)^2)
  r2_dec <- if (sstot > eps) 1 - ssres_dec/sstot else NA
  
  M <- r2_dec
  fit_y <- fit_dec
  direction <- "non-increasing"
  
  # residuals
  resid <- x - fit_y
  rdiff <- diff(resid)
  # raw residual jumpiness
  J_raw <- switch(residual_jump_fun,
                  MSSD = mean(rdiff^2),
                  MASD = mean(abs(rdiff)),
                  TotalVariation = sum(abs(rdiff)))
  
  # normalization factor
  norm_factor <- switch(normalize_by,
                        sd = if (sd(x) > eps) sd(x) else 1,
                        mad = if (mad(x, constant = 1) > eps) mad(x, constant = 1) else 1,
                        range = if (diff(range(x)) > eps) diff(range(x)) else 1,
                        none = 1)
  # if MSSD was used (squared), scale accordingly to get comparable units:
  if (residual_jump_fun == "MSSD") {
    # MSSD has squared units; normalize to variance for scale-free measure:
    # divide by var(x) (sd^2) if normalization is sd or mad or range chosen accordingly
    if (normalize_by == "sd") {
      J_norm <- J_raw / (norm_factor^2 + eps)
    } else if (normalize_by == "mad") {
      J_norm <- J_raw / (norm_factor^2 + eps)
    } else if (normalize_by == "range") {
      J_norm <- J_raw / ((norm_factor^2) + eps)
    } else {
      # if no normalization, keep raw MSSD but scale to something reasonable:
      J_norm <- J_raw
    }
  } else {
    # MASD or TotalVariation are in same units as x, so divide by norm_factor
    J_norm <- J_raw / (norm_factor + eps)
  }
  
  # Make J_norm non-negative; optionally cap or scale to [0,Inf). We leave as is.
  # Final metric
  MAJ <- (1 - M) * J_norm
  
  out <- list(
    n = n,
    isotonic_R2 = M,
    isotonic_direction = direction,
    residual_jump_raw = J_raw,
    residual_jump_norm = J_norm,
    MAJ = as.numeric(MAJ),
    details = list(residual_jump_fun = residual_jump_fun, normalize_by = normalize_by)
  )
  class(out) <- "monotone_adjusted_jumpiness"
  return(out)
}

print.monotone_adjusted_jumpiness <- function(x, digits = 4, ...) {
  cat("Monotone-Adjusted Jumpiness (MAJ)\n")
  cat("---------------------------------\n")
  cat("n:", x$n, "\n")
  cat("Isotonic R2 (best monotone fit):", formatC(x$isotonic_R2, digits = digits),
      " direction:", x$isotonic_direction, "\n")
  cat("Residual jump (raw):", formatC(x$residual_jump_raw, digits = digits),
      " normalized:", formatC(x$residual_jump_norm, digits = digits), "\n")
  cat("MAJ = (1 - R2) * normalized_residual_jumpiness :", formatC(x$MAJ, digits = digits), "\n")
  cat("Details:", paste(names(x$details), x$details, sep="="), "\n")
}














library(ggplot2)
library(patchwork)
library(energy)
# Idea: Make a bunch of cases of possible cumulative MAs
  # Try a bunch of metrics and see if one is what we want and makes sense
  # Let's do 7 studies
  # True mean of 0, true ts of 0.1, vi = 0.05 -> I^2 66%
set.seed(40)
studies <- c(1, 2, 3, 4, 5, 6, 7)

########################################################
# Case 1: Decreasing over time
########################################################
ts_1 <- c(0.4, 0.3, 0.25, 0.2, 0.14, 0.12, 0.1)
mean_1 <- rnorm(7, mean=0, sd=0.1)
pi_1 <- 2.571 * sqrt(ts_1 + 1/100)
lb_1 <- mean_1 - pi_1
ub_1 <- mean_1 + pi_1

case_1 <- as.data.frame(cbind(studies, mean_1, lb_1, ub_1))

case_1_MA <- ggplot(case_1, aes(x=studies, y=mean_1)) +
  geom_pointrange(aes(ymin=lb_1, ymax=ub_1), color="black", size=0.1)+
  geom_hline(yintercept=0) + 
  labs(x="Number of Studies Included",
       y="Theta") +
  scale_y_continuous(limits = c(-1.8, 1.8)) +
  theme_minimal()

case_1_ts <- ggplot(case_1, aes(x=studies, y=ts_1)) +
  geom_line(color="#07bbc1") + 
  geom_hline(yintercept=0.1) + 
  labs(x="Number of Studies Included", 
       y="TS") +
  scale_y_continuous(limits = c(-0.4, 0.4)) +
  theme_minimal()


########################################################
# Case 2: Increasing over time
########################################################
ts_2 <- c(0.1, 0.15, 0.16, 0.18, 0.2, 0.25, 0.27)
mean_2 <- rnorm(7, mean=0, sd=0.1)
pi_2 <- 2.571 * sqrt(ts_2 + 1/100)
lb_2 <- mean_2 - pi_2
ub_2 <- mean_2 + pi_2

case_2 <- as.data.frame(cbind(studies, mean_2, lb_2, ub_2))

case_2_MA <- ggplot(case_2, aes(x=studies, y=mean_2)) +
  geom_pointrange(aes(ymin=lb_2, ymax=ub_2), color="black", size=0.1)+
  geom_hline(yintercept=0) + 
  labs(x="Number of Studies Included",
       y="Theta") +
  scale_y_continuous(limits = c(-1.8, 1.8)) +
  theme_minimal()

case_2_ts <- ggplot(case_2, aes(x=studies, y=ts_2)) +
  geom_line(color="#07bbc1") + 
  geom_hline(yintercept=0.1) + 
  labs(x="Number of Studies Included", 
       y="TS") +
  scale_y_continuous(limits = c(-0.4, 0.4)) +
  theme_minimal()

########################################################
# Case 3: Oscillating big over time
########################################################
ts_3 <- c(0.1, 0.15, 0.04, 0.2, 0.07, 0.19, 0.1)
mean_3 <- rnorm(7, mean=0, sd=0.1)
pi_3 <- 2.571 * sqrt(ts_3 + 1/100)
lb_3 <- mean_3 - pi_3
ub_3 <- mean_3 + pi_3

case_3 <- as.data.frame(cbind(studies, mean_3, lb_3, ub_3))

case_3_MA <- ggplot(case_3, aes(x=studies, y=mean_3)) +
  geom_pointrange(aes(ymin=lb_3, ymax=ub_3), color="black", size=0.1)+
  geom_hline(yintercept=0) + 
  labs(x="Number of Studies Included",
       y="Theta") +
  scale_y_continuous(limits = c(-1.8, 1.8)) +
  theme_minimal()

case_3_ts <- ggplot(case_3, aes(x=studies, y=ts_3)) +
  geom_line(color="#07bbc1") + 
  geom_hline(yintercept=0.1) + 
  labs(x="Number of Studies Included", 
       y="TS") +
  scale_y_continuous(limits = c(-0.4, 0.4)) +
  theme_minimal()

########################################################
# Case 4: Oscillating little over time
########################################################
ts_4 <- c(0.1, 0.14, 0.06, 0.11, 0.08, 0.12, 0.1)
mean_4 <- rnorm(7, mean=0, sd=0.1)
pi_4 <- 2.571 * sqrt(ts_1 + 1/100)
lb_4 <- mean_4 - pi_4
ub_4 <- mean_4 + pi_4

case_4 <- as.data.frame(cbind(studies, mean_4, lb_4, ub_4))

case_4_MA <- ggplot(case_4, aes(x=studies, y=mean_4)) +
  geom_pointrange(aes(ymin=lb_4, ymax=ub_4), color="black", size=0.1)+
  geom_hline(yintercept=0) + 
  labs(x="Number of Studies Included",
       y="Theta") +
  scale_y_continuous(limits = c(-1.8, 1.8)) +
  theme_minimal()

case_4_ts <- ggplot(case_4, aes(x=studies, y=ts_4)) +
  geom_line(color="#07bbc1") + 
  geom_hline(yintercept=0.1) + 
  labs(x="Number of Studies Included", 
       y="TS") +
  scale_y_continuous(limits = c(-0.4, 0.4)) +
  theme_minimal()

########################################################
# Case 5: Decreasing then increasing
########################################################
ts_5 <- c(0.4, 0.25, 0.13, 0.1, 0.14, 0.2, 0.33)
mean_5 <- rnorm(7, mean=0, sd=0.1)
pi_5 <- 2.571 * sqrt(ts_5 + 1/100)
lb_5 <- mean_1 - pi_5
ub_5 <- mean_1 + pi_5

case_5 <- as.data.frame(cbind(studies, mean_5, lb_5, ub_5))

case_5_MA <- ggplot(case_5, aes(x=studies, y=mean_5)) +
  geom_pointrange(aes(ymin=lb_5, ymax=ub_5), color="black", size=0.1)+
  geom_hline(yintercept=0) + 
  labs(x="Number of Studies Included",
       y="Theta") +
  scale_y_continuous(limits = c(-1.8, 1.8)) +
  theme_minimal()

case_5_ts <- ggplot(case_5, aes(x=studies, y=ts_5)) +
  geom_line(color="#07bbc1") + 
  geom_hline(yintercept=0.1) + 
  labs(x="Number of Studies Included", 
       y="TS") +
  scale_y_continuous(limits = c(-0.4, 0.4)) +
  theme_minimal()

########################################################
# Case 6: Increasing then decreasing
########################################################
ts_6 <- c(0.05, 0.1, 0.16, 0.27, 0.19, 0.15, 0.1)
mean_6 <- rnorm(7, mean=0, sd=0.1)
pi_6 <- 2.571 * sqrt(ts_6 + 1/100)
lb_6 <- mean_6 - pi_6
ub_6 <- mean_6 + pi_6

case_6 <- as.data.frame(cbind(studies, mean_6, lb_6, ub_6))

case_6_MA <- ggplot(case_6, aes(x=studies, y=mean_6)) +
  geom_pointrange(aes(ymin=lb_6, ymax=ub_6), color="black", size=0.1)+
  geom_hline(yintercept=0) + 
  labs(x="Number of Studies Included",
       y="Theta") +
  scale_y_continuous(limits = c(-1.8, 1.8)) +
  theme_minimal()

case_6_ts <- ggplot(case_6, aes(x=studies, y=ts_6)) +
  geom_line(color="#07bbc1") + 
  geom_hline(yintercept=0.1) + 
  labs(x="Number of Studies Included", 
       y="TS") +
  scale_y_continuous(limits = c(-0.4, 0.4)) +
  theme_minimal()

########################################################
# Examine plots
########################################################
case_1_MA / case_1_ts
case_2_MA / case_2_ts
case_3_MA / case_3_ts
case_4_MA / case_4_ts
case_5_MA / case_5_ts
case_6_MA / case_6_ts

########################################################
# Metrics
# Spearman correlation - Measures monotonicity
# Kendall's tau - measures monotonicity, good for small samples with noise
# Monotonicity ratio - Fraction of differences that share majority direction (ie inc or dec)
# Isotonic regression R^2 
# Distance correlation - measures dependence via pairwise euclidean distances, measures strength of changes
# MSE
# Variance of the sequence
# MSSD - Captures how jumpy the seauence is
# MAJ
########################################################
diffs_1 <- diff(ts_1)
diffs_2 <- diff(ts_2)
diffs_3 <- diff(ts_3)
diffs_4 <- diff(ts_4)
diffs_5 <- diff(ts_5)
diffs_6 <- diff(ts_6)
eps = .Machine$double.eps

# Spearman
spearman_1 <- cor(studies, ts_1, method="spearman")
spearman_2 <- cor(studies, ts_2, method="spearman")
spearman_3 <- cor(studies, ts_3, method="spearman")
spearman_4 <- cor(studies, ts_4, method="spearman")
spearman_5 <- cor(studies, ts_5, method="spearman")
spearman_6 <- cor(studies, ts_6, method="spearman")

# Kendall
kendall_1 <- cor(studies, ts_1, method="kendall")
kendall_2 <- cor(studies, ts_2, method="kendall")
kendall_3 <- cor(studies, ts_3, method="kendall")
kendall_4 <- cor(studies, ts_4, method="kendall")
kendall_5 <- cor(studies, ts_5, method="kendall")
kendall_6 <- cor(studies, ts_6, method="kendall")

# Monotonicity ratio
mono_ratio_inc <- mean(diffs_1 >= 0)
mono_ratio_dec <- mean(diffs_1 <= 0)
mono_ratio_1 <- max(mono_ratio_inc, mono_ratio_dec)

mono_ratio_inc <- mean(diffs_2 >= 0)
mono_ratio_dec <- mean(diffs_2 <= 0)
mono_ratio_2 <- max(mono_ratio_inc, mono_ratio_dec)

mono_ratio_inc <- mean(diffs_3 >= 0)
mono_ratio_dec <- mean(diffs_3 <= 0)
mono_ratio_3 <- max(mono_ratio_inc, mono_ratio_dec)

mono_ratio_inc <- mean(diffs_4 >= 0)
mono_ratio_dec <- mean(diffs_4 <= 0)
mono_ratio_4 <- max(mono_ratio_inc, mono_ratio_dec)

mono_ratio_inc <- mean(diffs_5 >= 0)
mono_ratio_dec <- mean(diffs_5 <= 0)
mono_ratio_5 <- max(mono_ratio_inc, mono_ratio_dec)

mono_ratio_inc <- mean(diffs_6 >= 0)
mono_ratio_dec <- mean(diffs_6 <= 0)
mono_ratio_6 <- max(mono_ratio_inc, mono_ratio_dec)

# Isotonic Regresstion R^2
y <- ts_1
yhat_inc <- stats::isoreg(seq_along(y), y)$yf
ssres_inc <- sum((y - yhat_inc)^2)
sstot <- sum((y - mean(y))^2)
r2_inc <- if (sstot > eps) 1 - ssres_inc / sstot else NA
yhat_dec <- -stats::isoreg(seq_along(-y), -y)$yf
ssres_dec <- sum((y - yhat_dec)^2)
r2_dec <- if (sstot > eps) 1 - ssres_dec / sstot else NA
isotonic_r2_1 <- if (!is.na(r2_inc) && !is.na(r2_dec)) max(r2_inc, r2_dec) else NA

y <- ts_2
yhat_inc <- stats::isoreg(seq_along(y), y)$yf
ssres_inc <- sum((y - yhat_inc)^2)
sstot <- sum((y - mean(y))^2)
r2_inc <- if (sstot > eps) 1 - ssres_inc / sstot else NA
yhat_dec <- -stats::isoreg(seq_along(-y), -y)$yf
ssres_dec <- sum((y - yhat_dec)^2)
r2_dec <- if (sstot > eps) 1 - ssres_dec / sstot else NA
isotonic_r2_2 <- if (!is.na(r2_inc) && !is.na(r2_dec)) max(r2_inc, r2_dec) else NA

y <- ts_3
yhat_inc <- stats::isoreg(seq_along(y), y)$yf
ssres_inc <- sum((y - yhat_inc)^2)
sstot <- sum((y - mean(y))^2)
r2_inc <- if (sstot > eps) 1 - ssres_inc / sstot else NA
yhat_dec <- -stats::isoreg(seq_along(-y), -y)$yf
ssres_dec <- sum((y - yhat_dec)^2)
r2_dec <- if (sstot > eps) 1 - ssres_dec / sstot else NA
isotonic_r2_3 <- if (!is.na(r2_inc) && !is.na(r2_dec)) max(r2_inc, r2_dec) else NA

y <- ts_4
yhat_inc <- stats::isoreg(seq_along(y), y)$yf
ssres_inc <- sum((y - yhat_inc)^2)
sstot <- sum((y - mean(y))^2)
r2_inc <- if (sstot > eps) 1 - ssres_inc / sstot else NA
yhat_dec <- -stats::isoreg(seq_along(-y), -y)$yf
ssres_dec <- sum((y - yhat_dec)^2)
r2_dec <- if (sstot > eps) 1 - ssres_dec / sstot else NA
isotonic_r2_4 <- if (!is.na(r2_inc) && !is.na(r2_dec)) max(r2_inc, r2_dec) else NA

y <- ts_5
yhat_inc <- stats::isoreg(seq_along(y), y)$yf
ssres_inc <- sum((y - yhat_inc)^2)
sstot <- sum((y - mean(y))^2)
r2_inc <- if (sstot > eps) 1 - ssres_inc / sstot else NA
yhat_dec <- -stats::isoreg(seq_along(-y), -y)$yf
ssres_dec <- sum((y - yhat_dec)^2)
r2_dec <- if (sstot > eps) 1 - ssres_dec / sstot else NA
isotonic_r2_5 <- if (!is.na(r2_inc) && !is.na(r2_dec)) max(r2_inc, r2_dec) else NA

y <- ts_6
yhat_inc <- stats::isoreg(seq_along(y), y)$yf
ssres_inc <- sum((y - yhat_inc)^2)
sstot <- sum((y - mean(y))^2)
r2_inc <- if (sstot > eps) 1 - ssres_inc / sstot else NA
yhat_dec <- -stats::isoreg(seq_along(-y), -y)$yf
ssres_dec <- sum((y - yhat_dec)^2)
r2_dec <- if (sstot > eps) 1 - ssres_dec / sstot else NA
isotonic_r2_6 <- if (!is.na(r2_inc) && !is.na(r2_dec)) max(r2_inc, r2_dec) else NA

# Distance correlation
dcor_1 <- dcor(studies, ts_1)
dcor_2 <- dcor(studies, ts_2)
dcor_3 <- dcor(studies, ts_3)
dcor_4 <- dcor(studies, ts_4)
dcor_5 <- dcor(studies, ts_5)
dcor_6 <- dcor(studies, ts_6)

# MSE
mse_1 <- mean((ts_1 - 0.1)^2)
mse_2 <- mean((ts_2 - 0.1)^2)
mse_3 <- mean((ts_3 - 0.1)^2)
mse_4 <- mean((ts_4 - 0.1)^2)
mse_5 <- mean((ts_5 - 0.1)^2)
mse_6 <- mean((ts_6 - 0.1)^2)

# Variance of sequence
var_1 <- var(ts_1)
var_2 <- var(ts_2)
var_3 <- var(ts_3)
var_4 <- var(ts_4)
var_5 <- var(ts_5)
var_6 <- var(ts_6)

# MSSD
mssd_1   <- mean(diffs_1^2)
mssd_2   <- mean(diffs_2^2)  
mssd_3   <- mean(diffs_3^2)  
mssd_4   <- mean(diffs_4^2)  
mssd_5   <- mean(diffs_5^2)  
mssd_6   <- mean(diffs_6^2)  

# MASD
masd_1   <- mean(abs(diffs_1))
masd_2   <- mean(abs(diffs_2))   
masd_3   <- mean(abs(diffs_3))   
masd_4   <- mean(abs(diffs_4))   
masd_5   <- mean(abs(diffs_5))   
masd_6   <- mean(abs(diffs_6))   

# MAJ
maj_1 = monotone_adjusted_jumpiness(ts_1, residual_jump_fun = "MSSD", normalize_by="sd")[6]
maj_2 = monotone_adjusted_jumpiness(ts_2, residual_jump_fun = "MSSD", normalize_by="sd")[6]
maj_3 = monotone_adjusted_jumpiness(ts_3, residual_jump_fun = "MSSD", normalize_by="sd")[6]
maj_4 = monotone_adjusted_jumpiness(ts_4, residual_jump_fun = "MSSD", normalize_by="sd")[6]
maj_5 = monotone_adjusted_jumpiness(ts_5, residual_jump_fun = "MSSD", normalize_by="sd")[6]
maj_6 = monotone_adjusted_jumpiness(ts_6, residual_jump_fun = "MSSD", normalize_by="sd")[6]

# Combine
spearman = c(spearman_1, spearman_2, spearman_3, spearman_4, spearman_5, spearman_6)
kendall = c(kendall_1, kendall_2, kendall_3, kendall_4, kendall_5, kendall_6)
mono_ratio = c(mono_ratio_1, mono_ratio_2, mono_ratio_3, mono_ratio_4, mono_ratio_5, mono_ratio_6)
isotonic_r2 = c(isotonic_r2_1, isotonic_r2_2, isotonic_r2_3, isotonic_r2_4, isotonic_r2_5, isotonic_r2_6)
dcor = c(dcor_1, dcor_2, dcor_3, dcor_4, dcor_5, dcor_6)
mse = c(mse_1, mse_2, mse_3, mse_4, mse_5, mse_6)
variance = c(var_1, var_2, var_3, var_4, var_5, var_6)
mssd = c(mssd_1, mssd_2, mssd_3, mssd_4, mssd_5, mssd_6)
masd = c(masd_1, masd_2, masd_3, masd_4, masd_5, masd_6)
maj = c(maj_1, maj_2, maj_3, maj_4, maj_5, maj_6)
cases <- c(1, 2, 3, 4, 5, 6)

metrics <- as.data.frame(cbind(cases, spearman, kendall, mono_ratio, isotonic_r2, dcor, mse, variance, mssd, masd, maj))
rownames(metrics) <- NULL
metrics <- as.data.frame(lapply(metrics, as.numeric))
metrics <- round(metrics, 3)
tbl <- knitr::kable(metrics)
#save_kable(tbl, "metrics_table.txt")



(case_1_MA + case_2_MA) / (case_3_MA + case_4_MA) / (case_5_MA + case_6_MA)

(case_1_ts + case_2_ts) / (case_3_ts + case_4_ts) / (case_5_ts + case_6_ts)





