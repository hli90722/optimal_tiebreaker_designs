library(tidyverse)

# Inverse efficiency criterion
inv_eff <- function(Ez, Exz, Ex2z, Ex2){
  M11 <- 1-Ez^2-Exz^2/Ex2
  M12 <- -1*Exz*Ez-Ex2z*Exz/Ex2
  M22 <- Ex2-Exz^2-Ex2z^2/Ex2
  return(M11/(M11*M22-M12^2))
}

####################################### Helper functions
# Numerically compute value of E(x^2z) in [lower,upper] minimizing inv_eff 
# given Ez, Exz, Ex2
optimize_x2z <- function(Ez, Exz, Ex2, lower, upper, inv_eff,
                         tol=.Machine$double.eps){
  if (near(lower, upper, tol)){
    optimum <- list(minimum=lower,
                    objective=inv_eff(Ez, Exz, lower, Ex2))
  } else {
    optimum <- optimize(function(x) inv_eff(Ez, Exz, x, Ex2), 
                        lower=lower, upper=upper, tol=tol)
  }
  return(optimum)
}

# Computes \widetilde{xz}_{\max}(z_tilde)
# x_cumsum is a vector of cumulative sums of x when x is sorted
get_xz_max <- function(x_cumsum, z_tilde) {
  n <- length(x_cumsum)
  n_treat <- n * (1 + z_tilde) / 2
  whole_n_treat <- floor(n_treat)
  frac_n_treat <- n_treat - whole_n_treat
  whole_cumsum <- -1 * x_cumsum[(n - whole_n_treat)]
  x_thresh <- ifelse(whole_n_treat == 0, 0, x_cumsum[(n - whole_n_treat)]) - 
              ifelse(whole_n_treat == (n-1), 0, x_cumsum[(n - whole_n_treat - 1)])
  return(2 * (x_thresh * frac_n_treat + whole_cumsum) / n)
}

# Finds l = l(t, eps) where t = x[t_ind]
# This is the unique l such that E_p(z) = z_tilde for p(x) = l when x < t,
# p(t) = eps, and p(x) = 1 when x > t.
get_l <- function(t_ind, n, eps, z_tilde) {
  if (t_ind == 1) {
    return(0)
  }
  n_treat <- n * (1 + z_tilde) / 2
  return((n_treat - (n - t_ind) - eps) / (t_ind - 1))
}

# Finds eps such that E_p(xz) = xz_tilde for p(x) = l when x < t,
# p(t) = eps, and p(x) = 1 when x > t, where t = x[t_ind]
# x_cumsum is a vector of cumulative sums of x when x is sorted
get_eps_max <- function(t_ind, x_cumsum, z_tilde, xz_tilde) {
  lo_sum <- ifelse(t_ind == 1, 0, x_cumsum[(t_ind-1)])
  mid_sum <- x_cumsum[t_ind] - lo_sum
  hi_sum <- -1 * x_cumsum[t_ind]
  n <- length(x_cumsum)
  n_treat <- (1 + z_tilde) / 2 * n
  num <- n * xz_tilde / 2 - hi_sum - 
    lo_sum * (n_treat - (n - t_ind)) / (t_ind - 1) 
  den <- mid_sum - lo_sum / (t_ind - 1)
  return(num / den)
}

# Finds u = u(t, eps) where t = x[t_ind] for n sorted running variable values x.
# This is the unique u such that E_p(z) = z_tilde for p(x) = 0 when x < t
# p(t) = eps, and p(x) = u when x > t.
get_u <- function(t_ind, n, eps, z_tilde) {
  if (t_ind == n) {
    return(1)
  }
  n_treat <- n * (1 + z_tilde) / 2
  return((n_treat - eps) / (n - t_ind))
}

# Finds eps such that E_p(xz) = xz_tilde for p(x) = 0 when x < t,
# p(t) = eps, and p(x) = u when x > t, where t = x[t_ind] and u 
# is from get_u
# x_cumsum is a vector of cumulative sums of x when x is sorted + centered
get_eps_min <- function(t_ind, x_cumsum, z_tilde, xz_tilde) {
  lo_sum <- ifelse(t_ind == 1, 0, x_cumsum[(t_ind-1)])
  mid_sum <- x_cumsum[t_ind] - lo_sum
  hi_sum <- -1 * x_cumsum[t_ind]
  n <- length(x_cumsum)
  n_treat <- (1 + z_tilde) / 2 * n
  num <- n * xz_tilde / 2 - n_treat * hi_sum / (n - t_ind)
  den <- mid_sum - hi_sum / (n - t_ind)
  return(num / den)
}

# Computes E(x^az) for the 2 level design p(x) = l when x < t, p(t) = eps, and
# p(x) = u when x > t, where t = x[t_ind]
# Assumes x is sorted with cumulative sums of x^a  equal to cumsum_moments
compute_moments_from_cumsum <- function(cumsum_moments, l, u, t_ind, eps) {
  n <- length(cumsum_moments)
  lo_sum <- ifelse(t_ind == 1, 0, cumsum_moments[(t_ind-1)])
  mid_sum <- cumsum_moments[t_ind] - lo_sum
  hi_sum <- cumsum_moments[n] - cumsum_moments[t_ind]
  return((2 * (l * lo_sum + eps * mid_sum + u * hi_sum) - cumsum_moments[n]) / n)
}

# Evaluates the randomized two level design equal to l below t,
# eps at t, and u above t at x
eval_rand_two_level_design <- function(x, l, u, t, eps) {
  if (l > eps || eps > u) {
    return("Must have l <= eps <= u")
  }
  if (l < 0) {
    stop("Must have l >= 0")
  }
  if (t > 1) { 
    stop("Must have t <= 1")}
  return(l * (x < t) + eps * (x == t) + u * (x > t))
}

# Get l and u for a two level design, given cumulative sums of centered x,
# t_ind and eps=u, to satisfy constraints z_tilde, xz_tilde
get_l_u_high <- function(x_cumsum, t_ind, z_tilde, xz_tilde) {
  n <- length(x_cumsum)
  lo_sum <- ifelse(t_ind == 1, 0, x_cumsum[(t_ind-1)])
  n_treat <- (1 + z_tilde) / 2 * n
  A <- matrix(c(t_ind-1, lo_sum, n-t_ind+1, -1*lo_sum), nrow=2, byrow=FALSE)
  b <- c(n_treat, n*xz_tilde/2)
  soln <- solve(A) %*% b
  l <- soln[1]
  u <- soln[2]
  return(list(l=soln[1], u=soln[2]))
}

# Get l and u for a two level design, given cumulative sums of centered x,
# t_ind and eps=l, to satisfy constraints z_tilde, xz_tilde
get_l_u_lo <- function(x_cumsum, t_ind, z_tilde, xz_tilde) {
  n <- length(x_cumsum)
  lo_sum <- x_cumsum[t_ind]
  n_treat <- (1 + z_tilde) / 2 * n
  A <- matrix(c(t_ind, lo_sum, n-t_ind, -1*lo_sum), nrow=2, byrow=FALSE)
  b <- c(n_treat, n*xz_tilde/2)
  soln <- solve(A) %*% b
  l <- soln[1]
  u <- soln[2]
  return(list(l=soln[1], u=soln[2]))
}

# Get l and u for a two level design, given cumulative sums of centered x,
# t_ind, and eps, to satisfy constraints z_tilde, xz_tilde
get_l_u_given_eps <- function(x_cumsum, t_ind, eps, z_tilde, xz_tilde) {
  n <- length(x_cumsum)
  lo_sum <- ifelse(t_ind == 1, 0, x_cumsum[(t_ind-1)])
  mid_sum <- x_cumsum[t_ind] - lo_sum
  hi_sum <- -1 * x_cumsum[t_ind]
  n_treat <- (1 + z_tilde) / 2 * n
  A <- matrix(c(t_ind-1, lo_sum, n-t_ind, hi_sum), nrow=2, byrow=FALSE)
  b <- c(n_treat-eps, n*xz_tilde/2-eps*mid_sum)
  soln <- solve(A) %*% b
  l <- soln[1]
  u <- soln[2]
  return(list(l=soln[1], u=soln[2]))
}

# Function for root finder
# Given t_ind and eps, first uses get_l_u and x_cumsum (cumulative sums of x)
# to compute l and u to satisfy constraints given by z_tilde, xz_tilde.
# Then computes E(x^2z) using x2_cumsum (cumulative sums of x^2)
# Assumes x has been centered
get_x2z <- function(t_ind, eps, z_tilde, xz_tilde, x_cumsum, x2_cumsum) {
  lu <- get_l_u_given_eps(x_cumsum, t_ind, eps, z_tilde, xz_tilde)
  l <- lu$l
  u <- lu$u
  return(compute_moments_from_cumsum(x2_cumsum, l, u, t_ind, eps))
}

# Computes E[x^a*z] for 3 level tiebreaker
xaz_3level <- function(x, a, min_x, max_x){
  n <- length(x)
  p <- rep(0, n)
  p[x > max_x] <- 1
  p[(x >= min_x) & (x <= max_x)] <- 0.5
  return(mean(x^a*(2*p-1)))
}

#################################### Main functions
# Computes parameters for p_max_dag (if type="max") or p_min_dag (if type="min")
# given treatment fraction treat_frac 
# and normalized short-term gain parameter delta
get_p_dag <- function(x, treat_frac, delta, type="max", check_sorted=FALSE) {
  if (treat_frac <= 0 || treat_frac >= 1) {
    stop("treat_frac needs to be strictly between 0 and 1")
  }
  if (delta < 0 || delta > 1) {
    stop("delta needs to be between 0 and 1, inclusive")
  }
  if (type != "max" && type != "min"){
    stop("type needs to be either 'max' or 'min'")
  }
  x_mean <- mean(x)
  x <- x - x_mean
  n <- length(x)
  # Sort x if unsorted
  if (check_sorted && any(diff(x) < 0)) {
    x <- sort(x)
  }
  x_cumsum <- cumsum(x)
  
  # Convert from delta to xz_tilde
  z_tilde <- 2 * treat_frac - 1
  xz_max <- get_xz_max(x_cumsum, z_tilde)
  xz_tilde <- delta * xz_max
  if (type == "max") {
    # Binary search
    lo <- 1
    hi <- n + 1
    while (hi > (lo + 1)) {
      mid <- floor((lo + hi) / 2)
      test <- compute_moments_from_cumsum(cumsum_moments = x_cumsum,
                                          l = get_l(t_ind = mid, 
                                                    n = n, 
                                                    eps = 1, 
                                                    z_tilde = z_tilde),
                                          u = 1,
                                          t_ind = mid, 
                                          eps = 1)
      if (test >= xz_tilde) {
        lo <- mid
      } else {
        hi <- mid
      }
    }
    t <- x[lo]
    eps <- get_eps_max(lo, x_cumsum, z_tilde, xz_tilde)
    l <- get_l(lo, n, eps, z_tilde)
    return(list(l=l, eps=eps, t_ind=lo))
  } else {
    lo <- 0
    hi <- n
    while (hi > (lo + 1)) {
      mid <- floor((lo + hi) / 2)
      test <- compute_moments_from_cumsum(cumsum_moments = x_cumsum,
                                          l = 0,
                                          u = get_u(t_ind = mid, 
                                                    n = n, 
                                                    eps = 0, 
                                                    z_tilde = z_tilde),
                                          t_ind = mid, 
                                          eps = 0)
      if (test >= xz_tilde) {
        hi <- mid
      } else {
        lo <- mid
      }
    }
    t <- x[hi]
    eps <- get_eps_min(hi, x_cumsum, z_tilde, xz_tilde)
    u <- get_u(hi, n, eps, z_tilde)
    return(list(u=u, eps=eps, t_ind=hi))
  }
}

# Compute an optimal two-level design
# Uses optimize_x2z to numerically compute x^2z*(z_tilde, xz_tilde; eff) 
# where eff = inv_eff^{-1}
# If x is not mean centered, this function will operate on the centered x,
# then return parameters in terms of the original (uncentered) x
# lam_tol means if the computed lambda is within that amount of 0 or 1,
# we just return either p_max or p_min
get_opt_2_level <- function(x, treat_frac, delta, inv_eff, 
                            tol=.Machine$double.eps,
                            lam_tol=1e-5,
                            check_sorted=FALSE) {
  if (length(x) < 3) {
    stop("Error: running variable x must have length at least 3")
  }
  x_mean <- mean(x)
  if (check_sorted && any(diff(x) < 0)){
    x <- sort(x)
  }
  x <- x - x_mean
  x_cumsum <- cumsum(x)
  z_tilde <- 2 * treat_frac - 1
  xz_max <- get_xz_max(x_cumsum, z_tilde)
  xz_tilde <- delta * xz_max
  p_min_dag <- get_p_dag(x, treat_frac, delta, "min")
  p_max_dag <- get_p_dag(x, treat_frac, delta, "max")
  x2_cumsum <- cumsum(x^2)
  # Get endpoints of I_F from p_min_dag and p_max_dag
  min_x2z <- compute_moments_from_cumsum(cumsum_moments = x2_cumsum, 
                             l = 0,
                             u = p_min_dag$u,
                             t_ind = p_min_dag$t_ind,
                             eps = p_min_dag$eps)
  max_x2z <- compute_moments_from_cumsum(cumsum_moments = x2_cumsum, 
                             l = p_max_dag$l,
                             u = 1,
                             t_ind = p_max_dag$t_ind,
                             eps = p_max_dag$eps)
  x2z_opt <- optimize_x2z(Ez = z_tilde, 
                          Exz = xz_tilde, 
                          Ex2 = mean(x^2), 
                          lower = min_x2z, 
                          upper = max_x2z,
                          inv_eff = inv_eff,
                          tol = tol)$minimum
  lam <- ifelse(near(min_x2z, max_x2z, tol=tol),
                0,
                (max_x2z - x2z_opt))
  l <- NA
  u <- NA
  t <- NA
  eps <- NA
  if (near(lam, 0, tol=lam_tol)){
    return(list(l=p_max_dag$l, 
                u=1, 
                t_ind=p_max_dag$t_ind, 
                eps=p_max_dag$eps,
                lam=lam))
  } else if (near(lam, 1, tol=lam_tol)){
    return(list(l=0, 
                u=p_min_dag$u, 
                t_ind=p_min_dag$t_ind, 
                eps=p_min_dag$eps,
                lam=lam))
  } else {
    # loop over t_ind
    for (t_ind in p_min_dag$t_ind:p_max_dag$t_ind) {
      lower <- ifelse(t_ind==p_max_dag$t_ind, p_max_dag$eps,
                      get_l_u_lo(x_cumsum, t_ind, z_tilde, xz_tilde)$l)
      upper <- ifelse(t_ind==p_min_dag$t_ind, p_min_dag$eps,
                      get_l_u_high(x_cumsum, t_ind, z_tilde, xz_tilde)$u)
      failed <- FALSE
      soln <- try(soln <- uniroot(f=function(eps) get_x2z(t_ind, 
                                          eps, 
                                          z_tilde, 
                                          xz_tilde, 
                                          x_cumsum, 
                                          x2_cumsum) - x2z_opt,
                                  lower=lower,
                                  upper=upper,
                                  tol=tol),
                  silent = TRUE)
      if (!(class(soln) == "try-error")) {
        break
      }
    }
    eps <- soln$root
    lu <- get_l_u_given_eps(x_cumsum, t_ind, eps, z_tilde, xz_tilde)
    return(list(l=lu$l,
                u=lu$u,
                t_ind=t_ind,
                eps=eps,
                lam=lam))
  }
}

########################################## Head start example
# Read in data, sort by running variable x
hs <- read_csv("head_start.csv") %>%
  drop_na(povrate60) %>%
  arrange(povrate60)
x <- hs$povrate60 - mean(hs$povrate60)
x_cumsum <- cumsum(x)
n <- length(x)
Ex2 <- mean(x^2)
n_treat <- 300
treat_frac <- n_treat / n
z_tilde <- 2 * treat_frac - 1
xz_max <- get_xz_max(x_cumsum, z_tilde)

# Histogram of x
par(mfrow=c(1,2), mar=c(5,5,4,2))
hist(x, 
     freq=FALSE, 
     ylab="", 
     breaks=30, 
     main="", 
     xlab="Centered poverty index",
     cex.lab=1.7, 
     cex.axis=1.7, 
     cex.main=1.7, 
     cex.sub=1.7)
abline(v=x[2505], lty="dashed")

# Compute 3 level tie-breaker designs
min_x <- sapply(0:min(n_treat-1, n-n_treat-1), function(i) x[n-n_treat-i])
max_x <- sapply(0:min(n_treat-1, n-n_treat-1), function(i) x[n-n_treat+1+i])
xz_tiebreaker_3 <- sapply(1:length(min_x),
                          function(i) xaz_3level(x,1,min_x[i],max_x[i]))
delta_tiebreaker_3 <- xz_tiebreaker_3 / xz_max
x2z_tiebreaker_3 <- sapply(1:length(min_x),
                           function(i) xaz_3level(x,2,min_x[i],max_x[i]))
inv_eff_tiebreaker_3 <- sapply(1:length(min_x),
                               function(i) inv_eff(Ez=z_tilde, 
                                                   Exz=xz_tiebreaker_3[i], 
                                                   Ex2z=x2z_tiebreaker_3[i], 
                                                   Ex2=Ex2))
# Compute optimal monotone designs
delta_list <- seq(min(delta_tiebreaker_3), 1, length.out=100)
n_delta <- length(delta_list)
p_max_dag_list <- vector("list")
p_min_dag_list <- vector("list")
p_opt_dag_list <- vector("list")
for (i in 1:n_delta) {
  delta <- delta_list[i]
  p_max_dag_list[[i]] <- get_p_dag(x, treat_frac, delta, type="max")
  p_min_dag_list[[i]] <- get_p_dag(x, treat_frac, delta, type="min")
  p_opt_dag_list[[i]] <- get_opt_2_level(x, treat_frac, delta, inv_eff)
}
x2z_opt <- rep(NA, n_delta)
opt_inv_eff <- rep(NA, n_delta)
for (i in 1:n_delta) {
  x2z_opt[i] <- compute_moments_from_cumsum(cumsum_moments=cumsum(x^2),
                                            l=p_opt_dag_list[[i]]$l,
                                            u=p_opt_dag_list[[i]]$u,
                                            t_ind=p_opt_dag_list[[i]]$t_ind,
                                            eps=p_opt_dag_list[[i]]$eps)
  opt_inv_eff[i] <- inv_eff(Ez=z_tilde, 
                            Exz=delta_list[i]*xz_max, 
                            Ex2z=x2z_opt[i], 
                            Ex2=Ex2)
}

# Plotting
plot(delta_tiebreaker_3, 
     inv_eff_tiebreaker_3, 
     xlim=c(0.8, 1),
     ylim=c(0, max(inv_eff_tiebreaker_3)), 
     xlab=expression(paste("Normalized short-term gain ", delta)),
     # ylab=expression(Eff^{-1}),
     ylab=expression(Var(\hat{\beta}_3)),
     type="l",
     col="red",
     cex.lab=1.7, 
     cex.axis=1.7, 
     cex.main=1.7, 
     cex.sub=1.7,
     lwd=2)
lines(delta_list, opt_inv_eff, lwd=2, col="blue")
legend(x="topleft",
       legend=c("Three level tie-breaker", 
                "Optimal monotone design"),
       lty="solid",
       lwd=3,
       cex=1.7,
       col=c("red", "blue"))
