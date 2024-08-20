function(d1,d2=NULL,winterval=10,colour1="blue",colour2=NULL,pair=F,thr=0.05,yval=1){
  # Adds stars to indicate significance on top of plots created with plot.cis. 
  # d1 - data - (for one participant) or participants (for group-level analysis) x timepoints
  # d2 - more data
  # winterval - interval between successive timepoints in ms
  # colour1 - colour for plotting stars to indicate that d1 is significantly above chance, or, if d2 is set, significantly greater than d2
  # colour2 - colour for plotting stars to indicate that d2 is significantly greater than d1
  # pair - paired t-test (default = false)
  # thr - threshold for statistical significance (N.B. false discovery rate correction also applied)
  # yval - height to plot stars on figure
  
  # set number of timepoints
  ntimes=dim(d1)[2]
  # initialise store of p values and colours
  p <- rep(NA, times = ntimes)
  pc <- rep(colour1, times = ntimes)
  
  # conduct t-test.  If there
  # are two, test differences between them.
  for(i in c(1:ntimes)){
    # if there is one dataset, test difference from 0.5
    if (is.null(d2)){
      # if there is no variance, we cannot conduct a t-test, so simply set p to 1
      if(var(d1[,i])==0) p[i] <- 1 else{
      p[i] <- t.test(d1[,i],mu=0.5,alternative='greater')$p.value}
      # if there are two datasets, test difference between them
    } else {
      if(var(d1[,i])==0 | var(d2[,i])==0) p[i] <- 1 else{
      p[i] <- t.test(d1[,i],d2[,i],paired=pair)$p.value
      if(mean(d1[,i] - d2[,i]) > 0) pc[i] <- colour1 else pc[i] <- colour2}
    }
  }
  # control false discovery rate
  p <- p.adjust(p,method="fdr")
  
  # if there are at least some significant points
  if(sum(as.numeric(p <= thr), na.rm = T) > 0){
    # draw stars
    # - at the points where the p-value is significant
    # - at the height specified by yval
    # cex - size
    # - in the colour specified by pc
    # text((c(0:(ntimes-1))*winterval + 1)[p <= thr], rep(yval, times = ntimes)[p <= thr], labels = "*", cex = 1.5, col = pc[p <= thr])
    # alternatively, draw filled-in circles with the same parameters
    points((c(0:(ntimes-1))*winterval + 1)[p <= thr], rep(yval, times = ntimes)[p <= thr], pch = 16, cex = 1, col = pc[p <= thr])
    
    }
}