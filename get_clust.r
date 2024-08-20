function (mdat, k = 3, dmeth = "euc", cmeth = "hc"){
  # performs clustering. 
  # mdat - local temporal generalisation matrix
  # k - final number of clusters
  # dmeth - method for computing distance matrix. Options "euc" (Euclidean,
  # default), "cor" (correlation), "cos" (cosine).
  # cmeth - clustering method. Options "hc" (agglomerative hierarchical,
  # default), "km" (k-means)
  ## ENSURE THAT library(lsa) IS RUN BEFORE THIS FUNCTION IS CALLED.
  
  
  # transform the data depending on the option dmeth:
  # - "euc" - calculate the euclidean distance between each row of the matrix,
  # and make it squareform
  # - "cor" - calculate correlation distance between pairs of rows
  # - "cor" - calculate cosine distance between pairs of rows
  # t is used because cor and cos operate over columns by default and we need
  # then to operate over rows (i.e. train times).
  dmat <- switch(dmeth, euc = as.matrix(dist(mdat)), cor = 1-cor(t(mdat)), cos = 1-cosine(t(mdat)))
  # remove nans
  dmat[is.na(dmat)] <- 0
  # cluster and cut the tree into the specified number of clusters
  # as.dist takes the upper triangle
  dclust <- switch(cmeth, hc = cutree(hclust(as.dist(dmat)), 
      k), km = kmeans(as.dist(dmat), k)$cluster)
  # return 
  dclust
}
