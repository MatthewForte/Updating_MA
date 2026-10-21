library(energy)



###################################################
# MAJ Function
###################################################

# Monotone-Adjusted Jumpiness (MAJ)
monotone_adjusted_jumpiness <- function(x, 
                                        residual_jump_fun = c("MSSD","MASD","TotalVariation"),
                                        normalize_by = c("sd","mad","range","none"),
                                        eps = .Machine$double.eps) {
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

###################################################
# Get all metrics function
###################################################

get_all_metrics <- function(sequence) {
  
  diffs <- diff(sequence)
  
  # Spearman Correlation
  spearman <- cor(1:length(sequence), sequence, method="spearman")
  
  # Kendall's Tau
  kendall <- cor(1:length(sequence), sequence, method="kendall")
  
  # Monotonicity Ratio
  mono_ratio_inc <- mean(diffs >= 0)
  mono_ratio_dec <- mean(diffs <= 0)
  mono_ratio <- max(mono_ratio_inc, mono_ratio_dec)
  
  # Isotonic R2 (decreasing only)
  iso_dec <- isoreg(seq_along(sequence), -sequence)
  isotonic_r2 <- 1 - sum((sequence + iso_dec$yf)^2) / sum((sequence - mean(sequence))^2)
  
  # Distance Correlation
  dcor <- dcor(1:length(sequence), sequence)
  
  # Sample variance
  var <- var(sequence)
  
  # MSSD
  mssd <- mean(diffs^2)
  
  # MASD
  masd <- mean(abs(diffs))
  
  # MAJ
  maj_raw <- monotone_adjusted_jumpiness(sequence, residual_jump_fun="MSSD", normalize_by="none")
  maj_norm <- monotone_adjusted_jumpiness(sequence, residual_jump_fun="MSSD", normalize_by="sd")

  # Return  
  result <- list(
    spearman = spearman,
    kendall = kendall,
    mono_ratio = mono_ratio,
    isotonic_r2 = isotonic_r2,
    dcor = dcor,
    var = var,
    mssd = mssd,
    masd = masd,
    maj_raw = maj_raw[6],
    maj_norm = maj_norm[6]
  )
  return(result)
}









