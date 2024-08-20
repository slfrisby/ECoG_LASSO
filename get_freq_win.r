function(d, e = 1, ne = ne, w = w, ntpts = ntpts, nf = nf){
  # get_win loops over electrodes anyway, so here we just need to extract info at every frequency for a single electrode at a single timepoint.
  
  # make an index describing where the required data is located in the big data matrix
  tpindex <- ((e - 1) * ntpts) + w 
  freqindex <- ntpts*ne*(nf-1) 
  freqindex <- freqindex + tpindex
  # extract the columns from the data at each frequency specified by freqindex
  o <- d[,freqindex[1]]
  for (i in c(2:length(nf))) o <- cbind(o, d[,freqindex[i]])
  return(o)
}