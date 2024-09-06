# Detailed figures on the climate application

# get functions
source("functions.R")
source("paper_plotting.R")

# set colours (DW)
col1 = "darkmagenta"
col2 = "black"

col_disp1 = "deepskyblue2"
col_disp2 = "deepskyblue4"
col_shift1 = "orange1"
col_shift1_sat = "orange3"
col_shift2 = "orangered1"
col_shift2_sat = "orangered3"


# read in data and generate functions to evaluate densities and quantiles
dat1 <- read.csv("data_figures_climate/HadGEM2-ES_2x2.csv")
dF_model1 <- approxfun(density(dat1$model), yleft = 0, yright = 0)
dF_truth1 <- approxfun(density(dat1$truth), yleft = 0, yright = 0)
qF_model1 <- function(x) quantile(dat1$model, probs = x, type = 1)
qF_truth1 <- function(x) quantile(dat1$truth, probs = x, type = 1)


dat2 <- read.csv("data_figures_climate/HadGEM2-ES_3x3.csv")
dF_model2 <- approxfun(density(dat2$model), yleft = 0, yright = 0)
dF_truth2 <- approxfun(density(dat2$truth), yleft = 0, yright = 0)
qF_model2 <- function(x) quantile(dat2$model, probs = x, type = 1)
qF_truth2 <- function(x) quantile(dat2$truth, probs = x, type = 1)


# First figure: illustrating that different grid cells have biases into different directions
pdf("figure_climate2a.pdf", width = 4.5, height = 3)
ylim <- c(20, 50)
plot.densities.decomp(dF_model1, dF_truth1, qF_model1, qF_truth1, ylim = ylim,
                      use_legend = TRUE, lab_F = "M", lab_G = "E")
text(0.5, 48, "Example grid cell 1")
legend("bottomright", pch = c("M", "E"), legend = c("model", "empirical"),
       col = c(col1, col2), bty = "n")
dev.off()

pdf("figure_climate2b.pdf", width = 4.5, height = 3)
plot.densities.decomp(dF_model2, dF_truth2, qF_model2, qF_truth2, ylim = ylim,
                      use_legend = TRUE, lab_F = "M", lab_G = "E")
text(0.5, 48, "Example grid cell 2")
legend("bottomright", pch = c("M", "E"), legend = c("model", "empirical"),
       col = c(col1, col2), bty = "n")
dev.off()


# read in data for second example and generate functions
dat3 <- read.csv("data_figures_climate/disp_HadGEM2-ES_12x10.csv")
dF_model3 <- approxfun(density(dat3$model), yleft = 0, yright = 0)
dF_truth3 <- approxfun(density(dat3$truth), yleft = 0, yright = 0)
qF_model3 <- function(x) quantile(dat3$model, probs = x, type = 1)
qF_truth3 <- function(x) quantile(dat3$truth, probs = x, type = 1)

wd3 <- wd_decomp(qF_model3, qF_truth3)
# cd3 <- cd_decomp(qF_model3, qF_truth3) # too slow

# compute CD decomposition using qwcd_approx
alphas <- seq(0.005, 0.995, 0.005)
F = quantile(dat3$model, probs = alphas)
G = quantile(dat3$truth, probs = alphas)
cd3 <-  qwcd_approx(quantiles.F = F, quantiles.G = G,
            alphas = alphas, betas = alphas)

# Second figure: illustrating difference between AVM and CD
pdf("figure_climate3.pdf", width = 9/1.3, height = 3/1.3)
# par(las = 1)
layout(matrix(c(1, 2, 3, 4), nrow = 1), widths = c(1, 2.5, 2.5, 1))

ylim <- c(20, 50)
par(mar = c(2.2,2.2,0,0),mgp = c(1.2,0.4,0))
plot.densities(dF_model3,dF_truth3,colF = col1,colG = col2,ylim = ylim,ylab = "quantile",xaxt = "n")

par(mar = c(2.2,0,0,0),mgp = c(1.2,0.4,0))
plot.avm_decomp(qF_model3, qF_truth3, colF = col1, colG = col2, 
                ylim = ylim, yaxt = "n", lab_F = "M", lab_G = "E",
                use_legend = TRUE, legend_loc = "bottomleft")
text(0.5, 50, "Illustration of AVM")


plot.CD_decomp(qF_model3, qF_truth3, beta = 0.5, colF = col1, colG = col2, 
               ylim = ylim, yaxt = "n", lab_F = "M", lab_G = "E",
               use_legend = TRUE)
text(0.5, 50, expression(Illustration~of~CD~(beta==0.5)))

par(mar = c(2.2,0,0,2.5),mgp = c(1.2,0.4,0))
plot(NULL, xlim = c(0, 1), ylim = c(0, 3.5), axes = FALSE, 
     xlab = "")
mtext("divergence", side = 4, cex = 0.65, line = 1.5)
box()
axis(4)
text(0.3, 3.5, "AVM")
lines(c(0.3, 0.3), c(0, wd3[1]), col = col_shift1, lwd = 3)
lines(c(0.3, 0.3), c(wd3[1], sum(wd3)), col = col_disp1, lwd = 3)

text(0.7, 3.5, "CD")
lines(c(0.7, 0.7), c(0, cd3[2]), col = col_shift1, lwd = 3)
lines(c(0.7, 0.7), c(cd3[2], cd3[1]), col = col_disp1, lwd = 3)
dev.off()
