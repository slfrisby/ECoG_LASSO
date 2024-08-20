function (x=seq(0,1650,by=10), ymat=ymat, wwidth = 50, winterval = 10, rflag = F,colour) {
  # Plots AUC over time.
  # x - beginning of each window in ms
  # y - local temporal generalisation matrix
  # wwidth - width of window in ms
  # winterval - interval between successive timepoints in ms
  # rflag - Flag indicating whether r-squared of model should be plotted
  ############
  
  # library(flux) #Functions for plotting AUC
  #library(segmented) #Functions for fitting piecewise linear regression models
  
  # since we are interested in the area between the cluster timecourse and chance, subtract 0.5 from the LTG
  ymat <- ymat - 0.5
  y <- rep(0,times=)
  for (i in c(1:length(x))){
    y[i] <- auc(x,ymat[i,])
  }
  
  # fit model and get number of breakpoints
  tmp <- selgmented(lm(y~x),type="bic")
  nbp <- tmp$selection.psi$npsi
  
  nno <- 1650/wwidth # number of non-overlapping time windows
  xno <- seq(0,1650,by=wwidth) #Subset of non-overlapping windows beginning with the first
  yno <- y[match(xno,x)] #Same for y
  
  m <- segmented(lm(yno ~ xno), npsi = nbp) #Fit model with specified number of breakpoints
  
  #Plot model
  plot(m, ylim = c(-50,300), xlab = "Fit time (ms)", ylab = "Area under the timecourse", col = colour,
       lwd = 2, cex.lab= 1.5, cex.axis = 1.5)
  points(x, y, pch = 16, col = hsv(0, 0, 0.6, 0.5)) #Add all data points in gray
  points(xno, yno, pch = 16, cex = 1.2) #Add those used to fit model in black
  
  rsq = summary(m)$adj.r.square #Model r-square
  if (rflag) #Add to plot if flag set
    text(0, 300, labels = bquote(r^2 == .(rsq)), adj = 0, cex = 1.5)
  
  # add line at 0
  abline(h = 0, lty = 2)
  
}