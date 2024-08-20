function (mdat, clusters, clusterbysub, colour1, colour2, winterval=10){
  # Cluster classifiers and plots accuracy of each cluster over time.
  # mdat - local temporal generalisation matrix
  # clusters - vector indicating which classifiers belong to which cluster
  # (generated with get.clust)
  # clusterbysub - array of timecourses averaged by cluster (cluster x time x subjects)
  # winterval - interval between successive timepoints in ms
  ## ENSURE THAT library(lsa) IS RUN BEFORE THIS FUNCTION IS CALLED.
  
  # set colours
  cols <- colorRampPalette(c(colour1,colour2))
  # alternatively comment out this line and change "cols" to "rainbow" below
  
  # get number of classifiers and number of test times
  nclassif <- dim(mdat)[1]
  ntimes <- dim(mdat)[2]
  pdat <- matrix(0, max(clusters), ntimes)
  p <- matrix(0, max(clusters), ntimes)
  # plot
  # type "n" - do not plot data
  plot(1, 0.5, type = "n", ylim = c(0.4, 1), xlim = c(1,ntimes*winterval), 
       xlab = "Test time", ylab = "Accuracy", cex.lab = 1.5, cex.axis = 1.5)
  # plot dashed line for chance performance
  abline(h = 0.5, lty = 2)
  # for each cluster
  for (i in c(1:max(clusters))) {
    # if the cluster contains more than one classifier
    if (sum(clusters == i) > 1){
      # calculate the mean timecourse of accuracy for all classifiers in that cluster
      pdat[i,] <- colMeans(mdat[clusters == i, 1:ntimes])
      # else, if there is only one classifier in the cluster, store its timecourse 
    } else {pdat[i,] <- mdat[clusters == i, 1:ntimes]}
    # plot dashed lines to show cluster timecourses
    lines(c(0:(ntimes-1))*winterval, pdat[i,], lty = 5,  col = cols(max(clusters))[i])
    # test whether performance of that cluster is above chance and store p-value
    for (j in c(1:ntimes)){
      p[i,j] <- t.test(clusterbysub[i,j,],mu=0.5,alternative='greater')$p.value
    }
  }
  # control false discovery rate
  tmp <- c(p)
  tmp <- p.adjust(p,method="fdr")
  p <- matrix(tmp,nrow=max(clusters),ncol=ntimes)
  
  # plot solid lines where cluster performance is above chance following false
  # discovery rate correction (where it isn't, store NA to prevent plotting)
  pdat[p > 0.05] = NA
  for (i in c(1:max(clusters))){
    lines(c(0:(ntimes-1))*winterval, pdat[i,], lwd = 2,  col = cols(max(clusters))[i])
  }

  # plot dots to show which train times belong to each cluster. 
  # pch = 16 - plot points with filled-in circles
  points(c(0:(ntimes-1))*winterval, 0.9 + clusters/(max(clusters) * 10), pch = 16, col = cols(max(clusters))[clusters])
  
}