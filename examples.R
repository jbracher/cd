source("functions.R")

# Some distributions (distribution and quantile functions)
# Std. normal
pF.std.norm = pnorm
qF.std.norm = qnorm

# Other normal distributions
mean.shift = 2
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

# TODO: Other distributions?


# Approximated (quantile-weighted) Cramér distances
# Two distributions, equal, nice (symmetric, equidistant) quantile levels 
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


