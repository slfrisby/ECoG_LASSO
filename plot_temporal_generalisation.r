function (allpowerltg,ntpts= 166,colour1,colour2){
  # Plot temporal generalisation.
  
  # allpowerltg - matrix of train time x test time x subject
  # ntpts - number of timepoints
  # colour1 - colour to use to plot classifiers that perform above chance, but that perform significantly worse than the best classifier at that timepoint
  # colour2 - colour to use to plot classifiers that perform just as well as the best classifier at that timepoint
  
  # nclusters <- dim(tmp)[1]
  threshp <- matrix(0, ntpts, ntpts) # threshold p-value - is overall decoding significant at the time window?
  p <- matrix(0, ntpts, ntpts) #Output matrix containing p-values
  o <- matrix(0, ntpts, ntpts) #output matrix containing plotting data

  # check whether on-diagonal accuracy is reliably above threshold.
  # t.test defaults to checking whether the mean is different to zero
  for(i in c(1:ntpts)){
    for (j in c(1:ntpts)) {
      threshp[i,j] <- t.test(allpowerltg[i, j, ] - 0.5, alternative = "greater")$p.value
    }
  }
  # Control false discovery rate
  threshp <- matrix(p.adjust(threshp, method="fdr"), ntpts, ntpts)
  
  # Find best-performing classifier at each test timepoint. (Average over subjects first.)
  bestclassifier <- apply(apply(allpowerltg, c(1,2), mean), 2, which.max)
  
  # Compute contrast of all clusters with best-performing cluster at all test timepoints
   # For each test time (columns)
  for (j in c(1:ntpts)) {
    # Compute uncorrected p-value for contrast with best-performing classifier at each test timepoint
    # find the best classifier
    best <- bestclassifier[j]
    #For each classifier (rows)
    for (i in c(1:ntpts)) { 
      #If not the best
      if (i != best) { 
        # conduct a paired t-test to see whether accuracy for the best 
        # cluster at that timepoint is significantly greater than 
        # accuracy for cluster j at that timepoint
        p[i,j]<- t.test(allpowerltg[best, j, ], allpowerltg[i, j, ], paired = T, 
                        alternative = "greater")$p.value
        #If it is the best
      } else{ 
        p[i,j] <- 1.0
      }
    }
  }
  
  # Control false discovery rate
  p <- matrix(p.adjust(p, method="fdr"), ntpts, ntpts)
  
  #Fill in o with plotting data
  for(i in c(1:ntpts)){
    for (j in c(1:ntpts)) {
      # if decoding is not significant at all, plot grey
      if(threshp[i,j] > 0.05){
        o[i,j] <- 0
      } else { 
        # if the comparison p > 0.05, i.e. there is no significant difference between this
        # cluster and the best-performing cluster at the 5% significance level,
        # plot colour 2
        if (p[i,j] > 0.05)
          o[i,j] <- 2
        # if the comparison p < 0.05, i.e. there is a significant difference between this 
        # cluster and the best-performing cluster at the 0.5% significance
        # level, plot colour 1
        if (p[i,j] < 0.05)
          o[i,j] <- 1
      }
    }
  }
  
  # get vector of colours
  pcols <- c('#7D7F7C', colour1, colour2)
  
  # plot heatmap
  # xaxt = "n" - don't plot the x-axis
  # yaxt = "n" - don't plot the y-axis
  # image(t(o), xaxt = "n", yaxt = "n", col = pcols)
  image(seq(0,1650,by=10),seq(0,1650,by=10),o[,ntpts:1], xaxt = "n", yaxt = "n", col = pcols, xlab="Test time (ms) →",ylab="← Train time (ms)",cex.lab = 1.5)
  
  # draw a black border
  box()
  # axes
  axis(1, at = seq(0,1650,by=500),labels= c(0,500,1000,1500),cex.axis = 1.5)
  axis(2, at = seq(1650,0,by=-500),labels= c(0,500,1000,1500),cex.axis = 1.5)
  
  # don't print outputs
  invisible(list(o, p))
  
}