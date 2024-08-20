function (d, tp, w = 1, wsz = 50, ntpts = 166, nf = c(1:60)){
  # d - data (matrix of trials x spatiofrequency features (electrode 1 at 
  # timepoint 1 and frequency 1, electrode 1 at timepoint 2 and frequency 1, ...
  # electrode 1 at timepoint 34 and frequency 1, electrode 2 at timepoint 1 and 
  # frequency 1, ... electrode n at timepoint 34 and frequency 1, electrode 1 at 
  # timepoint 1 and frequency 2, ...)
  # tp - data type (frequencies "f", voltages "v")
  # wsz - window size (relevant only for voltage decoding)
  # nf = number of frequencies
	o <- NA
	if(tp=="f"){
	  ne <- dim(d)[2]/(ntpts*60) # 60 frequencies present in data
	  # for electrode 1 at timepoint 1, get matrix of items x frequencies
	  o <- get.freq.win(d, e = 1, ne = ne, w = w, ntpts = ntpts, nf = nf)
	  # for all other electrodes, get matrix of items x frequencies. Concatenate
	  # columns to produce one matrix containing data from all electrodes and all 
	  # frequencies at a single timepoint
	  if (ne>1){
  	  for (i in c(2:ne)){
  	    tmp <- get.freq.win(d, e = i, ne = ne, w = w, ntpts = ntpts, nf = nf)
  	    o <- cbind(o,tmp)}
	  }
		
		} else if(tp=="v"){
		    ne <- dim(d)[2]/4000 # 4000 timepoints present in the data (including baseline - stimulus onset is timepoint 1001)
		    # for electrode 1, get matrix of items x timepoints. The wsz/2 th timepoint is the centre of the wavelet used for frequencies 
		    # (e.g. for window 1, the wsz/2 th timepoint is stimulus onset)
		    o <- d[,c((1001+w-(0.5*wsz)):(1000+w+(0.5*wsz)))]
		    # for all other electrodes, get matrix of items x timepoints. Concatenate
		    # columns to produce one matrix of all electrodes and all timepoints in 
		    # a single window
		    if (ne > 1){
  		    for (i in c(2:ne)){
  		      tmp <- d[,c((1001+w-(0.5*wsz)+4000*(i-1)):(1000+w+(0.5*wsz)+4000*(i-1)))]
  		      o <- cbind(o,tmp)}
		    }
		    
		} else message("Unrecognized data type")
		o
}
