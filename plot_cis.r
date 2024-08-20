function(d, y.lim=c(0.4,1), newflag = T, winterval = 10, colour="blue"){
  # Compute mean and 95% confidence interval for each column of d
  # and plot as a ribbon plot.
  # d - results (folds x timepoints)
  # y.lim - y-axis limits 
  # newflag - indicates whether or not a new ribbon should be overlaid on the same plot
  # winterval - interval between successive timepoints in ms
  # colour - colour to plot
 
  ntimes <- dim(d)[2]
  # create results to plot - top row is mean, second row is lower confidence interval,
  # third row is upper confidence interval
	o <- matrix(0,3,ntimes)
	for(i in c(1:ntimes)){
		#If accuracy is constant, set mean, and set upper and lower CIs to same value
		if(var(d[,i])==0) o[,i] <- mean(d[,i]) else{	
			#Otherwise compute mean and CI
			t <- t.test(d[,i])
			o[1,i] <- t$estimate
			o[2:3,i] <- t$conf.int[1:2]
			}
		}
	
	# If newflag is set, initialise plot
	if(newflag){
	  # since no data is being plotted, the first 2 fields contain arbitrary values. 
	  # type "n" - do not plot data
		plot(1,0.5, type = "n", ylim=y.lim, xlim = c(1,ntimes*winterval), ylab="Accuracy", xlab="Time (ms)", cex.lab=1.5, cex.axis = 1.5)
	  # plot a line at 0.5 (chance)
	  # lty = 2 - plot dashed line
		abline(h = 0.5, lty = 2)
		}
	
	# Plot mean accuracy 
	x <- c(0:(ntimes-1))*winterval
	# lwd = 2 - line width 2
	lines(x, o[1,], type = "l", lwd=2, col = colour)
	# alternative version with points. 
	# pch = 16 - plot points with filled-in circles
	#lines(x, o[1,], type = "o", pch=16, lwd=2, col = colour)
	
	# Get pale colour for ribbon
	tmp <- tmp <- adjustcolor(colour,0.2)
	# Plot 95% confidence intervals
	polygon(c(x, x[ntimes:1]), c(o[2,], o[3,ntimes:1]), col=tmp, border=NA)
	invisible(o)
}