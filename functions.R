# Some functions to be used in our studies

############################################################
# Approximations of the Cramér distance
# Riemann sum approximation based on cumulative distribution functions
cd_CDF = function(pF,pG,step_x = 0.0001,lower = -100,upper = 100){
  # Do some safety checks? E.g., pF(lower),pG(lower) < eps, pF(upper),pG(upper) > 1 - eps ?
  grid_x = seq(from = lower, to = upper, by = step_x)
  return(sum(step_x*sum((pF(grid_x) - pG(grid_x))^2)))  
}

# Riemann sum approximation based on quantile functions
# Implement decomposition?
cd_quant = function(qF,qG, step_p = 0.0001){
  # Safety checks?
  alphas = betas = seq(from = step_p, to = 1 - step_p, by = step_p)
  gammas = c(0,alphas,1)
  K = length(alphas)
  
  quantiles.F = qF(alphas)
  quantiles.G = qG((gammas[-(K+2)] + gammas[-1])/2)
  
  integrand = function(i,j){
    return((gammas[j+1] - gammas[j])*(sign(alphas[i] - (gammas[j] + gammas[j+1])/2) != 
                                        sign(quantiles.F[i] - quantiles.G[j]))*
             abs(quantiles.F[i] - quantiles.G[j]))
  }
  return(2*sum(t(outer(1:K,1:(K+1),integrand)))/K)
}

############################################################
# Approximations of the quantile-weighted Cramér distance
# Approximation based on the full quantile function
# Implement decomposition?
qwcd = function(qF,qG,alphas,step_p = 0.0001,switch.args = FALSE){
  # Assumes uniform weights. Implement weighting?
  # Implement decomposition?
  if(switch.args){
    qG.prior = qG
    qG = qF
    qF = qG.prior
  }
  
  gammas = seq(from = 0, to = 1, by = step_p)
  K = length(alphas)
  L = length(gammas)
  
  quantiles.F = qF(alphas)
  quantiles.G = qG((gammas[-L] + gammas[-1])/2)
  
  integrand = function(i,j){
    return((gammas[j+1] - gammas[j])*
             (sign(alphas[i] - (gammas[j] + gammas[j+1])/2) != 
                sign(quantiles.F[i] - quantiles.G[j]))*
             abs(quantiles.F[i] - quantiles.G[j]))
  }
  return(2*sum(t(outer(1:K,1:(L-1),integrand)))/K)
}

# Approximations based on a finite number of known quantiles
# Simple approximation formula for the (quantile-weighted) Cramér distance
# (only special cases of two distribution setting)
# TODO: Implement further safety checks + decomposition
qwcd_simple = function(quantiles.F = NULL, quantiles.G = NULL, qF = NULL, qG = NULL,alphas = NULL,betas = NULL){
  if(is.null(alphas)) stop("Quantile levels (alphas) need to be specified!")
  K = length(alphas)
  if(is.null(betas)){
    betas = alphas
    print("Quantile levels of G (betas) set to match quantile levels of F (alphas).")
  }
  
  if(is.null(quantiles.F) + is.null(qF) != 1) 
    stop("The first distribution needs to be specified either through a vector of quantiles 
         at the given levels (alphas) or through its quantile function (qF). Do not specify both!")
  if(is.null(quantiles.G) + is.null(qG) + is.null(y) != 2)
    stop("The second distribution needs to be specified either through a vector of quantiles 
         at the given levels (alphas) or through its quantile function (qG). Do not specify both!")
  
  quantiles.F = qF(alphas)
  quantiles.G = qG(betas)
  
  integrand = function(i,j) (sign(alphas[i] - alphas[j]) != sign(quantiles.F[i] - quantiles.G[j]))*
    abs(quantiles.F[i] - quantiles.G[j])/(K+1)
  return(2*sum(t(outer(1:K,1:K,integrand)))/K)
}

# General approximation formula for the quantile-weighted Cramér distance + decomposition
# (empirical setting and two distribution setting)
# TODO: Implement a method to symmetrize quantile levels, e.g., by linear interpolation
qwcd_approx = function(quantiles.F = NULL, quantiles.G = NULL, qF = NULL, qG = NULL, y = NULL,
                       alphas = NULL, betas = NULL, weights = NULL, return_decomp = TRUE){
  if(is.null(alphas)) stop("Quantile levels (alphas) need to be specified!")
  K = length(alphas)
  # Check required symmetries: Quantile levels (alphas) should be symmetric around 0.5 
  # for the decomposition to make sense.
  if(any(round(alphas + alphas[K:1],10) != 1))
    warning("Quantile levels do not bound central prediction intervals.")
  
  if(is.null(quantiles.F) + is.null(qF) != 1) 
    stop("The first distribution needs to be specified either through a vector of quantiles 
         at the given levels (alphas) or through its quantile function (qF). Do not specify both!")
  if(is.null(quantiles.G) + is.null(qG) + is.null(y) != 2)
    stop("The second distribution needs to be specified either through a vector of quantiles 
         at the given levels (betas or, if betas are not specified, alphas), through its quantile 
         function (qF) or through a sample (y). Do not use more than one way of specifying G!")
  
  if(is.null(quantiles.F)){
    quantiles.F = qF(alphas)
    print("Quantiles of F computed from quantile function.")
  }
  else if(length(quantiles.F) != length(alphas)) 
    stop("Number of quantiles of F (quantiles.F) does not match number of quantile levels
         (alphas).")
  
  if(is.null(y)){
    if(is.null(betas)){
      betas = alphas
      print("Quantile levels of G (betas) set to match quantile levels of F (alphas).")
    }
    else{
      M = length(betas)
      if(any(round(betas + betas[M:1],10) != 1)) warning("Quantile levels (betas) not symmetric.")
    }
    if(is.null(quantiles.G)){
      quantiles.G = qG(betas)
      print("Quantiles of G computed from quantile function.")
    }
    else if(length(quantiles.G) != length(betas)) 
      stop("Number of quantiles of G (quantiles.G) does not match number of quantile levels
           (betas/alphas).")
    qGhat = approxfun(x = c(0,betas,1), y = c(min(quantiles.F,quantiles.G), quantiles.G,
                                              max(quantiles.F,quantiles.G)))
  }
  else{
    qGhat = qG = function(p) quantile(y,p = p,type = 1)
    M = length(y)
    betas = seq(1/M,1,1/M)
    print("Empirical quantile function used as second distribution (G).")
  }
  
  if(is.null(weights)){
    weights = rep(1/K,K)
    print("No weights specified. Use uniform weights.")
  }
  else{ # Check symmetry.
    if(any(weights != weights[K:1])) 
      warning("Weights are not symmetric. The decomposition only uses the first half 
            of the weight vector to weight the central prediction intervals and 
            therefore does not match the qwCD.")
  }
  
  gammas = sort(unique(round(c(0,alphas,betas,1),digits = 10))) 
  # Rounding avoids multiples that differ only by an epsilon
  L = length(gammas) # This is L+1 in (2)!
  
  # Compute components
  integrand_comps = function(i,j){
    weight_gamma = gammas[j+1] - gammas[j]
    lF = quantiles.F[i]
    uF = quantiles.F[K+1 - i]
    lG = qGhat((gammas[j+1] + gammas[j])/2)
    uG = qGhat((gammas[L+1 - (j+1)] + gammas[L+1 - j])/2)
    
    return(ifelse(alphas[i] == 0.5 || round((gammas[j+1] + gammas[j])/2,10) == 0.5, 1, 2)*
             # correction factor for the median times factor 2
             weights[i]*weight_gamma*                   # Assumes symmetric weights for simplicity...
             c(SFG = pmax(0,pmin(lF - lG, uF - uG) + pmax(0,lF - uG)),
               SGF = pmax(0,pmin(lG - lF, uG - uF) + pmax(0,lG - uF)),
               DFG = ifelse(gammas[j] + gammas[j+1] <= 2*alphas[i],
                            pmax(0, (uF - lF) - (uG - lG)), 0),
               DGF = ifelse(2*alphas[i] <= gammas[j] + gammas[j+1],
                            pmax(0, (uG - lG) - (uF - lF)), 0)))
  }
  
  comps = rowSums(apply(expand.grid(1:ceiling(K/2),1:ceiling((L-1)/2)), 1,
                        function(x) integrand_comps(x[1],x[2])))
  
  # Compute approximation directly
  # Can be skipped to speed-up computation. Replace with:
  # qwcd = sum(comps)
  
  integrand_qwcd = function(i,j){
    weight_gamma = (gammas[j+1] - gammas[j])
    quantile.Ghat = qGhat((gammas[j+1] + gammas[j])/2)
    
    return(weights[i]*weight_gamma*abs(quantiles.F[i] - quantile.Ghat)*
             (sign(alphas[i] - (gammas[j+1] + gammas[j])/2) != sign(quantiles.F[i] - quantile.Ghat)))
  }
  qwcd = 2*sum(t(outer(1:K,1:(L-1),integrand_qwcd)))
  
  # Check: Approximated qwCD is equal sum of components
  if(round(qwcd,10) != round(sum(comps),10)) 
    warning(paste0("Sum of components does not match approximate qwCD. 
                   Difference (qwCD - Sum) = ",signif((qwcd - sum(comps)))))
  
  if(return_decomp) return(c(qwcd = qwcd, SFG = comps[['SFG']], SGF = comps[['SGF']],
                             DFG = comps[['DFG']], DGF = comps[['DGF']]))
  else return(qwcd)
}