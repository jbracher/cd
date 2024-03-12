source("functions.R")
source("plotting_routines.R")

options(scipen = 9999)

# Quantile function of a mixture of uniform distributions
qF.pw.unif = function(p,bounds = c(-3,-1,0,1,3),weights = NULL){
  n = length(bounds) - 1
  if(!all(bounds[-1] >= bounds[-(n+1)])) stop("Misspecified bounds detected.")
  if(is.null(weights)) weights = rep(1/n,n)
  if(n != length(weights)) stop("Wrong number of weights!")
  else weights = weights/sum(weights)
  cumweights = c(0,cumsum(weights))
  n_p = length(p)
  bin = rowSums(matrix(rep(cumweights,each = n_p),nrow = n_p) <= p)
  x1 = bounds[bin]
  x2 = bounds[bin+1]
  return(x1 + (x2-x1)*(p - cumweights[bin])/weights[bin])
}

# Corresponding density function (if there are no point masses)
dF.pw.unif = function(x,bounds = c(-3,-1,0,1,3),weights = NULL){
  n = length(bounds) - 1
  if(!all(bounds[-1] > bounds[-(n+1)])) stop("Point mass or misspecified bounds detected.")
  if(is.null(weights)) weights = rep(1/n,n)
  if(n != length(weights)) stop("Wrong number of weights!")
  else weights = weights/sum(weights)
  cumweights = c(0,cumsum(weights))
  n_x = length(x)
  bin = rowSums(matrix(rep(bounds,each = n_x),nrow = n_x) <= x)
  x1 = bounds[ifelse(bin == 0,1,bin)]
  x2 = bounds[bin+1]
  return(ifelse(bin == 0 | bin == n+1,0,weights[ifelse(bin == 0,1,bin)]/(x2-x1)))
}


# Examples
# Figure 3 / Example 3.1
qF = function(p) qF.pw.unif(p,bounds = c(-2,0,2))
qG = function(p) qF.pw.unif(p,bounds = c(-2,0,1))

# pdf("figs/diffShift_centralInts.pdf", width = 3,height = 3)
# par(mar = c(2,2,0,0) + 0.2,mgp = c(1.2,0.4,0))
plot_avm_decomp(qF,qG)
# dev.off()

wd_decomp(qF,qG)/wd(qF,qG) # 1-WD: 100% dispersion
pwd_decomp(qF,qG)/pwd(qF,qG) # squared 2-WD
cd_decomp(qF,qG)/sum(cd_decomp(qF,qG)) # CD: 50% shift


# Figure 4 / Example 3.4
qH = function(p) qF.pw.unif(p,c(-5,0,1,5),c(3,2,1))
qF = function(p) qH(p)
qG = function(p) 2*qH(p)

# pdf("figs/CD-shift_loc-scale.pdf", width = 3,height = 3)
# par(mar = c(2,2,0,0) + 0.2,mgp = c(1.2,0.4,0))
plot_avm_decomp(qF,qG)
# dev.off()

cd_decomp(qF,qG)/sum(cd_decomp(qF,qG))
wd_decomp(qF,qG)/wd(qF,qG)
pwd_decomp(qF,qG,2)/pwd(qF,qG,2)


# Figure 5 / Example 3.6
qF = function(p) qF.pw.unif(p,bounds = c(-1,-1,0,0),weights = c(1,0,3))
qG = function(p,loc = 0.5,scale = 1.8) qF(p)*scale + loc

# pdf("figs/pWD_loc-scale.pdf", width = 3,height = 3)
# par(mar = c(2,2,0,0) + 0.2,mgp = c(1.2,0.4,0))
plot_avm_decomp(qF,qG)
# dev.off()

wd_decomp(qF,qG)/wd(qF,qG)
pwd_decomp(qF,qG,2)/pwd(qF,qG,2)
pwd_decomp(qF,qG,3)/pwd(qF,qG,3)


# Figure 6 / Example 3.8
qF = function(p) qF.pw.unif(p,bounds = 10*c(-0.1,0.1,0.2),weights = c(1,1))
qG = function(p) qF.pw.unif(p,bounds = 10*c(-0.05,0,0.2),weights = c(1,1))

# pdf("figs/pWD_unimodal.pdf", width = 3,height = 3)
# par(mar = c(2,2,0,0) + 0.2,mgp = c(1.2,0.4,0))
plot_avm_decomp(qF,qG)
# dev.off()

wd_decomp(qF,qG)/wd(qF,qG)
pwd_decomp(qF,qG,2)/pwd(qF,qG,2)
pwd_decomp(qF,qG,3)/pwd(qF,qG,3)


# Figure 7 / Example 3.9
qF = function(p) qF.pw.unif(p,bounds = c(-1,-1,0.1,1.2,1.2),weights = c(3,2,2,3))
qG = function(p) qF.pw.unif(p,bounds = c(-1,-1,0,1,1),weights = c(3,2,2,3))

# pdf("figs/CD-WD_sym.pdf", width = 3,height = 3)
# par(mar = c(2,2,0,0) + 0.2,mgp = c(1.2,0.4,0))
plot_avm_decomp(qF,qG)

wd_decomp(qF,qG)/wd(qF,qG) # 1-WD: 20% shift
pwd_decomp(qF,qG)/pwd(qF,qG) # squared 2-WD
cd_decomp(qF,qG)/sum(cd_decomp(qF,qG)) # CD: 8.6% shift


# Figure 8 / Example 3.11
qF = function(p) qF.pw.unif(p,bounds = c(-3,-2,-1,0,1,2),weights = c(1,0,1,1,1))
qG = function(p) qF.pw.unif(p,bounds = c(-2,-1,0,1,2,3),weights = c(1,1,1,0,1))

# pdf("figs/dispOrder_centralInts.pdf", width = 3,height = 3)
# par(mar = c(2,2,0,0) + 0.2,mgp = c(1.2,0.4,0))
plot_avm_decomp(qF,qG)
# dev.off()

cd_decomp(qF,qG)/sum(cd_decomp(qF,qG))


# Figure 9 / Example 3.15
qF = function(p) qF.pw.unif(p,bounds = c(-4,-1,-0.1,0,3),weights = c(2,0.1,1,3.1))
qG = function(p) qF.pw.unif(p,bounds = c(-1,0,1))
qH = function(p) qF.pw.unif(p,bounds = c(-3,0,0.1,1,4),weights = c(3.1,1,0.1,2))

# pdf("figs/intransitiveCD_centralInts.pdf", width = 3,height = 3)
# par(mar = c(2,2,0,0) + 0.2,mgp = c(1.2,0.4,0))
plot_avm_decomp(qF,qG)
plot_avm_decomp(qG,qH,colF = "red",colG = "green")
plot_avm_decomp(qF,qH,colG = "green")
# dev.off()

cd_decomp(qF,qG)/sum(cd_decomp(qF,qG))
cd_decomp(qG,qH)/sum(cd_decomp(qG,qH)) 
cd_decomp(qF,qH)/sum(cd_decomp(qF,qH)) 


# Figure 10 / Example 3.20
qF = function(p) qF.pw.unif(p,bounds = c(-4,-1.8,0,0.5,3,4),weights = c(4,6,5,1,4))
qG = function(p) qF.pw.unif(p,bounds = c(-4,0,4))
qH = function(p) qF.pw.unif(p,bounds = c(-4,-3,-0.5,0,1.8,4),weights = c(4,1,5,6,4))

# pdf("figs/intransitive_centralInts.pdf", width = 3,height = 3)
# par(mar = c(2,2,0,0) + 0.2,mgp = c(1.2,0.4,0))
plot_avm_decomp(qF,qG)
plot_avm_decomp(qG,qH,colF = "red",colG = "green")
plot_avm_decomp(qF,qH,colG = "green")
# dev.off()

wd_decomp(qF,qG)/wd(qF,qG)
wd_decomp(qG,qH)/wd(qG,qH)
wd_decomp(qF,qH)/wd(qF,qH)

cd_decomp(qF,qG)/sum(cd_decomp(qF,qG))
cd_decomp(qG,qH)/sum(cd_decomp(qG,qH)) 
cd_decomp(qF,qH)/sum(cd_decomp(qF,qH)) 



