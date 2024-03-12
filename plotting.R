# Some plotting routines to visualize the decompositions

# Visualization of the AVM decomposition
plot_avm_decomp = function(qF,qG,colF = "blue",colG = "red",col_dispG = "indianred2",col_dispF = "orange",col_shiftG = "hotpink1",col_shiftG_double = "deeppink3",col_shiftF = "violet",col_shiftF_double = "darkviolet"){
  avm_lines <- function(l1, u1, l2, u2, at, # interval ends and an x coordinate (at)
                        col_disp1 = rgb(1, 0.5, 0, 0.3), # colours
                        col_disp2 = rgb(0, 1, 0, 0.3),
                        col_shift1 = rgb(1, 0, 0, 0.3),
                        col_shift1_double = rgb(1, 0, 0, 0.6),
                        col_shift2 = rgb(0, 0, 1, 0.3),
                        col_shift2_double = rgb(0, 0, 1, 0.6), ...){
    
    # compute interval widths
    d1 <- u1 - l1
    d2 <- u2 - l2
    
    # initialize vectors where line segments will be stored
    disp1 <- disp2 <-
      shift1 <- shift1_double <-
      shift2 <- shift2_double <- NA
    
    # case where interval 1 is wider
    if(d1 >= d2){
      d <- d1 - d2
      
      # run through cases
      if(u1 < l2){
        disp1 <- c(l1, l1 + d)
        shift2 <- c(u2, l2, NA, u1, l1 + d)
        shift2_double <- c(l2, u1)
      }
      
      if(l1 <= l2 & l2 < u1 & u1 <= u2){
        disp1 <- c(l1, l1 + d)
        shift2 <- c(l1 + d, l2, NA, u1, u2)
      }
      
      if(l1 <= l2 & u2 <= u1){
        disp1 <- c(l1, l2, NA, u2, u1)
      }
      
      if(l2 <= l1 & l1 < u2 & u2 <= u1){
        disp1 <- c(u1, u1 - d)
        shift1 <- c(u1 - d, u2, NA, l1, l2)
      }
      
      if(u2 < l1){
        disp1 <- c(u1, u1 - d)
        shift1 <- c(u1 - d, l1, NA, u2, l2)
        shift1_double <- c(l1, u2)
      }
    }
    
    # case where interval 2 is wider
    if(d2 > d1){
      d <- d2 - d1
      
      # run through cases
      if(u2 < l1){
        disp2 <- c(l2, l2 + d)
        shift1 <- c(u1, l1, NA, u2, l2 + d)
        shift1_double <- c(l1, u2)
      }
      
      if(l2 <= l1 & l1 < u2 & u2 <= u1){
        disp2 <- c(l2, l2 + d)
        shift1 <- c(l2 + d, l1, NA, u2, u1)
      }
      
      if(l2 <= l1 & u1 <= u2){
        disp2 <- c(l2, l1, NA, u1, u2)
      }
      
      if(l1 <= l2 & l2 < u1 & u1 <= u2){
        disp2 <- c(u2, u2 - d)
        shift2 <- c(u2 - d, u1, NA, l2, l1)
      }
      
      if(u1 < l2){
        disp2 <- c(u2, u2 - d)
        shift2 <- c(u2 - d, l2, NA, u1, l1)
        shift2_double <- c(l2, u1)
      }
    }
    
    # plot lines
    lines(rep(at, length(disp1)), disp1, col = col_disp1, ...)
    lines(rep(at, length(disp2)), disp2, col = col_disp2, ...)
    lines(rep(at, length(shift1)), shift1, col = col_shift1, ...)
    lines(rep(at, length(shift1_double)), shift1_double, col = col_shift1_double, ...)
    lines(rep(at, length(shift2)), shift2, col = col_shift2, ...)
    lines(rep(at, length(shift2_double)), shift2_double, col = col_shift2_double, ...)
  }
  
  # function to plot illustration of area validation metric
  plot_avm <- function(l1, u1, l2, u2, alpha, # lower and upper interval ends, levels
                       col1 = "red", # colours
                       col2 = "blue",
                       col_disp1 = rgb(1, 0.5, 0, 0.3),
                       col_disp2 = rgb(0, 1, 0, 0.3),
                       col_shift1 = rgb(1, 0, 0, 0.3),
                       col_shift1_double = rgb(1, 0, 0, 0.6),
                       col_shift2 = rgb(0, 0, 1, 0.3),
                       col_shift2_double = rgb(0, 0, 1, 0.6),
                       legend_intervals = TRUE, legend_components = TRUE, ylim = NULL, ...){ # additional graphcal parameters
    # compute suitable ylim if not provided
    if(is.null(ylim)) ylim <- range(c(l1, l2, u1, u2))
    plot(NULL, type = "l", xlim = 0:1, ylim = ylim,
         xlab = "coverage", ylab = "interval limits")
    
    # run through values in alpha and plot lines
    for(i in 1:length(alpha)){
      avm_lines(l1[i], u1[i], l2[i], u2[i], at = alpha[i], lwd = 2,
                col_disp1 = col_disp1, col_disp2 = col_disp2,
                col_shift1 = col_shift1, col_shift1_double = col_shift1_double,
                col_shift2 = col_shift2, col_shift2_double = col_shift2_double)
    }
    
    # add lines for interval ends
    lines(alpha, l1, col = col1, lwd = 2,lty = 3)
    lines(alpha, u1, col = col1, lwd = 2, lty = 5)
    lines(alpha, l2, col = col2, lwd = 2, lty = 3)
    lines(alpha, u2, col = col2, lwd = 2, lty = 5)
    
    # add legends
    if(legend_intervals){
      legend("topleft", col = c(NA, col1, col1, NA, col2, col2),
             legend = c("Interval ends:", "upper F", "lower F",
                        NA, "upper G", "lower G"),
             lty = c(NA, 5, 3, NA, 5, 3), lwd = 2, ncol = 2, bty = "n")
    }
    if(legend_components){
      legend("bottomleft", col = c(NA, col_disp1, col_shift1, col_shift1_double, NA, col_disp2, col_shift2, col_shift2_double),
             legend = c("Components:", "F more dispersed", "F shifted up", "F shifted up (counted twice)",
                        NA, "G more dispersed", "G shifted up", "G shifted up (counted twice)"), pch = 15, bty = "n", ncol = 2)
    }
  }
  
  alpha <- 1:200/201
  
  # compute interval ends:
  l2 <- qF((1 - alpha)/2)
  u2 <- qF(1 - (1 - alpha)/2)
  l1 <- qG((1 - alpha)/2)
  u1 <- qG(1 - (1 - alpha)/2)
  
  # plot:
  plot_avm(l1, u1, l2, u2,col1 = colG,col2 = colF, alpha,col_disp1 = col_dispG,col_disp2 = col_dispF,col_shift1 = col_shiftG,col_shift1_double = col_shiftG_double,col_shift2 = col_shiftF,col_shift2_double = col_shiftF_double,legend_intervals = FALSE,legend_components = FALSE)
}

