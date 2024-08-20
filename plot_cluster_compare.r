function (tmp=clustacc.by.sj,aa = allaccs,ntpts= 166,colour1,colour2){
#Given array coding the mean decoding accuracy of each cluster
#at every time point for each subject, generates a heat plot
#indicating where the best-performing cluster is reliably better than
#other clusters.
#tmp: Array containing mean accuracy data for each subject
#thresh: Only compare clusters for time-points where mean accuracy exceeds this threshold
#allaccs: Array containing LTG matrices for each subject

  nclusters <- dim(tmp)[1]
  o <- matrix(0, nclusters, ntpts) #output matrix containing significance level
	p <- matrix(0, nclusters, ntpts) #Output matrix containing actual probabilities
  threshp <- rep(0,ntpts) # threshold p-value - is overall decoding significant at the time window?
  
  # check whether on-diagonal accuracy is reliably above threshold.
  # t.test defaults to checking whether the mean is different to zero
	for (i in c(1:ntpts)){
	  threshp[i] <- t.test(aa[i, i, ] - 0.5, alternative = "greater")$p.value
	}
  # control false discovery rate
  threshp <- p.adjust(threshp,method="fdr")

	# Find best-performing cluster at each timepoint. (Average over subjects first.)
	bestclust <- apply(apply(tmp, c(1,2), mean), 2, which.max)
	
	# Compute contrast of all clusters with best-performing cluster at all timepoints
	# For each window
  for (i in c(1:ntpts)) {
    # if on-diagonal accuracy is not above threshold, set values to NA
    if (threshp[i] > 0.05) {
 		  p[, i] <- NA
      # Otherwise compute uncorrected p-value for contrast with best-performing cluster for each cluster
      } else { 
        # find the best cluster
        best <- bestclust[i]
        #For each cluster
        for (j in c(1:nclusters)) { 
          #If not the best
          if (j != best) { 
            # conduct a paired t-test to see whether accuracy for the best 
            # cluster at that timepoint is significantly greater than 
            # accuracy for cluster j at that timepoint
            p[j,i]<- t.test(tmp[best, i, ], tmp[j, i, ], paired = T, 
                    alternative = "greater")$p.value
            #If it is the best
            } else{ 
				      p[j,i] <- 1.0
				    }
          }
       }
    }
		
	# Control false discovery rate
	p <- matrix(p.adjust(p, method="fdr"), nclusters, ntpts)

	#Fill in o with p-values, binned by significance band 
	for(i in c(1:ntpts)){
		for (j in c(1:nclusters)) {
			if(is.na(p[j,i])){
			  # if the value is NA, i.e. decoding is not significant at all, plot black
				o[j,i] <- 0
				} else { 
				  # if p > 0.005, i.e. there is no significant difference between this 
				  # cluster and the best-performing cluster at the 0.5% significance level
				  # (but there is at the 1% significance level), plot blue
  			  if (p[j,i] > 0.005) 
  				o[j, i] <- 2
  			  # if p > 0.01, i.e. there is no significant difference between this
  			  # cluster and the best-performing cluster at the 1% significance level
  			  # (but there is at the 5% significance level), plot green
  			  if (p[j,i] > 0.01)
  				o[j, i] <- 3
  			  # if p > 0.05, i.e. there is no significant difference between this
  			  # cluster and the best-performing cluster at the 5% significance level,
  			  # plot yellow
  			  if (p[j,i] > 0.05)
  				o[j, i] <- 4
  			  # if p <= 0.005, i.e. there is a significant difference between this 
  			  # cluster and the best-performing cluster at the 0.5% significance
  			  # level, plot purple
  			  if (p[j,i] <= 0.005)
  				o[j, i] <- 1
				}
			}
		}
	
	# get vector of colours
	#pcols <- viridis(4)
	colours <- colorRampPalette(c(colour1,colour2))
	pcols <- colours(4)
	# add black as the first colour
	pcols <- c('#000000FF', pcols)
	
	# plot heatmap
	# xaxt = "n" - don't plot the x-axis
	# yaxt = "n" - don't plot the y-axis
  # image(t(o), xaxt = "n", yaxt = "n", col = pcols)
	image(seq(0,1650,by=10),seq(0,1650,by=10),o[,ntpts:1], xaxt = "n", yaxt = "n", col = pcols, xlab="Test time (ms) →",ylab="← Train time (ms)",cex.lab = 1.5)
	
  # draw a black border
  box()
  # label axis with cluster number
  # side = 2 - write text down the left-hand side
  # las = 2 - write text perpendicular to the axis (i.e. the right way up)
  # line = 0.5 - specification of distance of text from plot
  # cex = font size
  # at = provides coordinates of each string of text
  # mtext(side = 2, las = 2, line = 0.5, cex = 1, at = c(0:(nclusters-1)/(nclusters-1)), 
      #text = c(1:nclusters))
  # alternatively, label axis with time
  axis(1, at = seq(0,1650,by=500),labels= c(0,500,1000,1500),cex.axis = 1.5)
  axis(2, at = seq(1650,0,by=-500),labels= c(0,500,1000,1500),cex.axis = 1.5)
  
  
  # don't print outputs
  invisible(list(o, p))
}
