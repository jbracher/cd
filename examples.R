source("functions.R")
source("check_formulas_cd.R")

# Some distributions (distribution and quantile functions)
# Std. normal
pF.std.norm = pnorm
qF.std.norm = qnorm

# Other normal distributions
mean.shift = 1
sd.disp = 2

# Shifted from std. normal
pF.shift.norm = function(q) pnorm(q,mean = mean.shift)
qF.shift.norm = function(p) qnorm(p,mean = mean.shift)

# More dispersed than std. normal
pF.disp.norm = function(q) pnorm(q,sd = sd.disp)
qF.disp.norm = function(p) qnorm(p,sd = sd.disp)

# Shift + Dispersion
pF.shift.disp.norm = function(q) pnorm(q,mean = mean.shift,sd = sd.disp)
qF.shift.disp.norm = function(p) qnorm(p,mean = mean.shift,sd = sd.disp)

# Uniform distribution
low = -1.6
up = 1.6
pF.unif = function(q) punif(q,min = low,max = up)
qF.unif = function(p) qunif(p,min = low,max = up)

# Mixtures of uniforn distributions
pF.pwlin1 = function(q) ifelse(q < 0,0,ifelse(q < 4,0.5*q/4,ifelse(q < 6, 0.5 + (q-4)/4, 1)))
pF.pwlin2 = function(q) ifelse(q < 1,0,ifelse(q < 3,(q-1)/4,ifelse(q < 7, 0.5 + (q-3)/8, 1)))

qF.pwlin1 = function(p) ifelse(p < 0.5, 8*p, 4 + (p-0.5)*4)
qF.pwlin2 = function(p) ifelse(p < 0.5, 1 + 4*p, 3 + (p-0.5)*8)


# Approximated (quantile-weighted) Cramér distances
# Two distributions, equal, nice (symmetric, equidistant) quantile levels
compare_values = function(pF,qF,pG,qG,alphas){
  round(rbind(c(cd_CDF(pF,pG),cd_decomp(qF,qG)),
        c(qwcd(qF,qG,alphas), qwcd_decomp(qF,qG,alphas)),
        c(qwcd(qG,qF,alphas), qwcd_decomp(qG,qF,alphas))[c(1,3,2,5,4)],
        unlist(approx_cd_v3(qF(alphas),qG(alphas))),
        qwcd_approx(qF = qF,qG = qG,alphas = alphas),
        qwcd_approx(qF = qG,qG = qF,alphas = alphas)[c(1,3,2,5,4)],
        cd_approx(qF(alphas),qG(alphas),alphas),
        cd_approx(qF(alphas),qG(alphas),alphas,add_levels = seq(0.01,0.99,0.01))),4)
}

alphas = seq(0.1,0.9,0.1)

compare_values(pF.disp.norm,qF.disp.norm,pF.std.norm,qF.std.norm,alphas)
compare_values(pF.shift.norm,qF.shift.norm,pF.std.norm,qF.std.norm,alphas)
compare_values(pF.shift.disp.norm,qF.shift.disp.norm,pF.std.norm,qF.std.norm,alphas)
compare_values(pF.unif,qF.unif,pF.std.norm,qF.std.norm,alphas)
compare_values(pF.pwlin1,qF.pwlin1,pF.pwlin2,qF.pwlin2,alphas)


alphas = seq(0.1,0.9,0.1)

pF = pF.std.norm
pG = pF.shift.norm
qF = qF.std.norm
qG = qF.shift.norm

# or
pF = pF.std.norm
pG = pF.shift.disp.norm
qF = qF.std.norm
qG = qF.shift.disp.norm

# or
pF = pF.shift.norm
pG = pF.disp.norm
qF = qF.shift.norm
qG = qF.disp.norm

c(cd_CDF(pF,pG),
  cd_quant(qF,qG),
  qwcd(qF,qG,alphas),
  qwcd(qF,qG,alphas,switch.args = TRUE),
  qwcd_simple(qF = qF,qG = qG,alphas),
  qwcd_approx(qF = qF,qG = qG,alphas = alphas))

# TODO: More examples
# Empirical setting, unequal quantile levels, etc.


