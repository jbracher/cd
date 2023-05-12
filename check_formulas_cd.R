# Version 1 of the approximation: summing over penalty terms; equation \eqref{eq:approx_cd}
# q_F: vector containing the (1:K)/(K + 1) quantiles of F
# q_G: vector containing the (1:K)/(K + 1) quantiles of G
approx_cd_v1 <- function(q_F, q_G){
  K <- length(q_F)
  grid_p <- (1:K)/(K + 1) # function assumes that the quantile levels are equally spaced
  
  contributions_F <- rep(NA, length(grid_p))
  for(i in seq_along(grid_p)){
    contributions_F[i] <- sum(
      (sign(grid_p[i] - grid_p) != sign(q_F[i] - q_G))*abs(q_F[i] - q_G)
    )
  }
  
  cd <- 2*sum(contributions_F)/K/(K + 1)
  return(cd)
}

# Version 2 of the approximation: summing over "almost squared" differences in CDF; equation \eqref{eq:rectangles}
# q_F: vector containing the (1:K)/(K + 1) quantiles of F
# q_G: vector containing the (1:K)/(K + 1) quantiles of G
approx_cd_v2 <- function(q_F, q_G){
  
  # compute quantile levels from length of provided quantile vectors:
  K <- length(q_F)
  if(length(q_G) != K) stop("q_F and q_G need to be of the same length")
  p <- (1:K)/(K + 1) # function assumes that the quantile levels are equally spaced
  
  # pool quantiles:
  q0 <- c(q_F, q_G)
  # vector of grouping variables, with 1 for values belonging to F, -1 for values 
  # belonging to G
  a0 <- c(rep(1, length(q_F)), rep(-1, length(q_G)))
  
  # re-order both vectors:
  q <- q0[order(q0)]
  a <- a0[order(q0)]
  # and compute "how many quantiles ahead" F or G is at a given segment:
  b <- abs(cumsum(a))
  
  # compute the lengths of segments defined by sorted quantiles:
  diffs_q <- c(diff(q), 0) # zero necessary for indexing below, but we could put 
  # anything (gets multiplied w zero)
  
  # and approximate CD
  cvm <- sum(diffs_q*b*(b + 1))/(K + 1)/(K)
  
  return(mean(cvm))
}

# implementation of the interval divergence:
# l_F: lower ends of prediction intervals for F
# u_F: upper ends of prediction intervals for F
# alpha_F: nominal coverage levels of intervals for F
# l_G, u_G, alpha_G: same for G
interval_comparison <- function(l_F, u_F, alpha_F, l_G, u_G, alpha_G){
  # computation of components via a distinction of cases:
  
  # if both PIs have same level
  if(alpha_F == alpha_G){
    F_disp <- max((u_F - l_F) - (u_G - l_G), 0)
    G_disp <- max((u_G - l_G) - (u_F - l_F), 0)
    F_larger <- max(max(u_F - u_G, 0) + max(l_F - l_G, 0) + max(l_F - u_G, 0) - F_disp - G_disp, 0)
    G_larger <- max(max(u_G - u_F, 0) + max(l_G - l_F, 0) + max(l_G - u_F, 0) - F_disp - G_disp, 0)
  }
  # if F has lower nominal coverage and "should" be nested in G:
  if(alpha_F < alpha_G){
    F_disp <- max((u_F - l_F) - (u_G - l_G), 0)
    G_disp <- 0
    F_larger <- max(max(u_F - u_G, 0) + max(l_F - u_G, 0) - F_disp, 0)
    G_larger <- max(max(l_G - l_F, 0) + max(l_G - u_F, 0) - F_disp, 0)
  }
  # if G has lower nominal coverage and "should" be nested in F:
  if(alpha_G < alpha_F){
    G_disp <- max((u_G - l_G) - (u_F - l_F), 0)
    F_disp <- 0
    G_larger <- max(max(u_G - u_F, 0) + max(l_G - u_F, 0) - G_disp, 0)
    F_larger <- max(max(l_F - l_G, 0) + max(l_F - u_G, 0) - G_disp, 0)
  }
  
  id <- F_larger + G_larger + F_disp + G_disp
  
  return(list(id = id,
              F_larger = F_larger,
              G_larger = G_larger,
              F_disp = F_disp,
              G_disp = G_disp))
}

# Version 3 of the approximation: via interval divergence; equation \eqref{eq:via_intervals}
# q_F: vector containing the (1:K)/(K + 1) quantiles of F
# q_G: vector containing the (1:K)/(K + 1) quantiles of G
approx_cd_v3 <- function(q_F, q_G){
  # compute quantile levels from length of provided quantile vectors:
  K <- length(q_F)
  if(length(q_G) != K) stop("q_F and q_G need to be of the same length")
  p <- (1:K)/(K + 1) # function assumes that the quantile levels are equally spaced
  coverages <- 2*(1:(K/2))/(K + 1)
  
  n_intervals <- K/2  
  # matrices to store differences and booleans indicating whether there is a discrepancy and
  # in which component it should be counted
  matrix_interval_comparisons <-
    matrix_F_larger <- matrix_G_larger <-
    matrix_F_disp <- matrix_G_disp <-
    matrix(NA, ncol = n_intervals, nrow = n_intervals,
           dimnames = list(paste("F", 1:n_intervals), paste("G", 1:n_intervals)))
  
  # fill these matrices:
  for(i in 1:n_intervals){
    for(j in 1:n_intervals){
      i_comp <- interval_comparison(l_F = q_F[i], u_F = q_F[K + 1 - i], alpha_F = p[K + 1 - i] - p[i],
                                     l_G = q_G[j], u_G = q_G[K + 1 - j], alpha_G = p[K + 1 - j] - p[j])
      matrix_interval_comparisons[i, j] <- i_comp$id
      matrix_F_larger[i, j] <- i_comp$F_larger
      matrix_G_larger[i, j] <- i_comp$G_larger
      matrix_F_disp[i, j] <- i_comp$F_disp
      matrix_G_disp[i, j] <- i_comp$G_disp
    }
  }
  
  cd <- 2*sum(matrix_interval_comparisons)/K/(K + 1)
  F_larger <- 2*sum(matrix_F_larger)/K/(K + 1)
  G_larger <- 2*sum(matrix_G_larger)/K/(K + 1)
  F_disp <- 2*sum(matrix_F_disp)/K/(K + 1)
  G_disp <- 2*sum(matrix_G_disp)/K/(K + 1)
  
  return(list(cd = cd, 
              F_larger = F_larger, G_larger = G_larger,
              F_disp = F_disp, G_disp = G_disp))
}


# independent implementation of the weighted interval score
weighted_interval_score_v1 <- function(q_F, observed){
  K <- length(q_F)
  levels <- (1:K)/(K + 1)
  2*mean(((observed < q_F) - levels)*(q_F - observed))
}

# independent implementation of the weighted interval score incl decomposition from previous project
# q_F: vector containing the (1:K)/(K + 1) quantiles of F
# observed: the observation y
# weights: optional weight for the different interval levels
# detailed: should a more detialed output be generated?
weighted_interval_score_v2 <- function(q_F, observed, weights = NULL, detailed = FALSE){
  
  K <- length(q_F)
  alpha <- 1 - (2*1:(K/2) - 1)/(K + 1)
  
  if(is.null(weights)){
    weights <- alpha/2
    weights[alpha == 1] <- weights[alpha == 1]/2 # weigh down interval corresponding to absolute error as the IS corresponds to 2 times the AE
  }
  
  if(length(weights) != length(alpha)){
    stop("weights and alpha need to be of the same length.")
  }
  
  alpha_half <- alpha/2
  one_m_alpha_half <- 1 - alpha/2
  
  # compute relevant quantiles:
  l <- rev(head(q_F, K/2))
  u <- tail(q_F, K/2)
  
  # compute widths of prediction intervals:
  widths_pi <- u - l
  
  # compute penalties:
  penalties_l <- 2/alpha*pmax(0, l - observed)
  penalties_u <- 2/alpha*pmax(0, observed - u)
  
  # compute interval scores at different levels:
  interval_scores <- widths_pi + penalties_l + penalties_u
  
  # name vectors
  names(l) <- names(u) <- names(widths_pi) <- names(penalties_l) <- names(penalties_u) <-
    names(interval_scores) <- names(alpha) <- names(alpha_half) <-
    names(one_m_alpha_half) <- paste0("alpha.", alpha)
  
  # compute combined score:
  numerator <- length(alpha) - 0.5*as.numeric(1 %in% alpha) # count ae only half if contained
  
  weighted_penalty_l <- sum(weights*penalties_l)/numerator
  weighted_penalty_u <- sum(weights*penalties_u)/numerator
  weighted_width_pi <- sum(weights*widths_pi)/numerator
  weighted_interval_score <- sum(weights*interval_scores)/numerator
  
  if(detailed){
    return(list(
      l = l, u = u, observed = observed,
      width_pi = widths_pi, penalty_l = penalties_l, penalty_u = penalties_u,
      interval_score = interval_scores, weights = weights,
      weighted_penalty_l = weighted_penalty_l,
      weighted_penalty_u = weighted_penalty_u,
      weighted_width_pi = weighted_width_pi,
      weighted_interval_score = weighted_interval_score
    ))
  }else{
    list(weighted_penalty_l = weighted_penalty_l,
         weighted_penalty_u = weighted_penalty_u,
         weighted_width_pi = weighted_width_pi,
         weighted_interval_score = weighted_interval_score)
  }
}

# Example:
K <- 1000
p <- (1:K)/(K + 1) # quantile levels - using many so approximation is good

# define parameters of two normals:
# F
mu_F <- 12
sigma_F <- 5
# G
mu_G <- 9
sigma_G <- 4

# compute quantiles:
q_F <- qnorm(p, mu_F, sigma_F) # quantiles of F
q_G <- qnorm(p, mu_G, sigma_G) # quantiles of G

# check all three implementations give exactly the same result:
approx_cd_v1(q_F, q_G)
approx_cd_v2(q_F, q_G)
approx_cd_v3(q_F, q_G)
# correct

# in addition compare to sampling-based approximation:
n_sim <- 10^6
X <- rnorm(n_sim, mean = mu_F, sd = sigma_F)
X_dash <- rnorm(n_sim, mean = mu_F, sd = sigma_F)

Y <- rnorm(n_sim, mean = mu_G, sd = sigma_G)
Y_dash <- rnorm(n_sim, mean = mu_G, sd = sigma_G)

mean(abs(X - Y)) - 0.5*(mean(abs(X - X_dash)) + mean(abs(Y - Y_dash)))

# compare to weighted interval score, which should be a special case:
observed <- 3
approx_cd_v1(q_F, rep(observed, K))
approx_cd_v3(q_F, rep(observed, K))

weighted_interval_score_v1(q_F, observed)
weighted_interval_score_v2(q_F, observed)
# everything checks out, incl decomposition


