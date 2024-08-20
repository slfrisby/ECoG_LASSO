function(d,ntpts=166){
  # make plots of time x frequency (for one electrode at a time). 
  # d - data (time x frequency)
  
  # library(fields)
  
  c <- matrix(0,dim(d)[1],dim(d)[2])
  
  # categorise coefficients - red for theta, orange for alpha, yellow for beta, light green for gamma, dark green for high gamma
  for (i in c(1:dim(d)[1])){
    for (j in c(1:dim(d)[2])){
      if (d[i,j] > 0 | d[i,j] < 0){ 
        if (j <= 11){ c[i,j] <- 1
        } else if (j <= 18){ c[i,j] <- 2
        } else if (j <= 31){ c[i,j] <- 3
        } else if (j <= 41){ c[i,j] <- 4
        } else {c[i,j] <- 5}
      }
    }
  }

  # plot 
  colours = c("#FFFFFF","#A50026","#F46D43","#FDCC3F","#66BD63","#006837")
  image(c(1:dim(d)[1]), c(1:dim(d)[2]), c, axes = F, xlim = c(1,ntpts), ylim = c(1,60), zlim=c(0,5),xlab = "Time (ms)", ylab = "Frequency (Hz)", col = colours, cex.lab = 1.5)
  box()
  axis(1,seq(1,ntpts,by=50),labels= c(0,500,1000,1500),cex.axis = 1.5)
  axis(2,seq(1,60,by=10),labels=c(0,10,20,30,40,50),cex.axis = 1.5)
  
  ### OLD ###
  # categorise coefficients - red for positive, blue for negative, green for neutral
  # for (i in c(1:dim(d)[1])){
  #   for (j in c(1:dim(d)[2])){
  #     if (d[i,j] > 0) { c[i,j] <- 0.5
  #     } else if (d[i,j] < 0) { c[i,j] <- -0.5
  #     } else {c[i,j] <- 0}
  #   }
  # }
  # 
  # # plot
  # image(c(1:dim(d)[1]), c(1:dim(d)[2]), c, axes = F, xlim = c(1,ntpts), ylim = c(1,60), zlim=c(-1,1),xlab = "Time (ms)", ylab = "Frequency (Hz)", col = rev(rainbow(3)), cex.lab = 1.5)
  # axis(1,seq(0,ntpts,by=5),labels=seq(0,1650,by=50), cex.axis = 1.5)
  # axis(2,seq(0,60,by=5), cex.axis = 1.5)
  # # plot white lines to divide bands. The squares indicating a coefficient are 
  # # CENTRED on the y tick they correspond to - the .5 ensures that the lines 
  # # divide the bands used for the main analysis.
  # abline(h = 11.5, col = "white")
  # abline(h = 18.5, col = "white")
  # abline(h = 31.5, col = "white")
  # abline(h = 41.5, col = "white")
}