function(d, tp, y, wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:60), return.models=T){
  # d - data (matrix of trials x spatiofrequency features (electrode 1 at 
  # time -1000 and frequency 1, electrode 1 at time -9999 and frequency 1, ...
  # electrode 1 at time 2999 and frequency 1, electrode 2 at timepoint -1000 and 
  # frequency 1, ... electrode n at time 2999 and frequency 1, electrode 1 at 
  # time -1000 and frequency 2, ...) or spatiotemporal features (electrode 1 at 
  # time -1000, electrode 1 at time -999, ... electrode 1 at time 2999, 
  # electrode 2 at time -1000) 
  # tp - data type (frequencies "f", voltages "v")
  # y - labels vector. For living/nonliving data, averaged across repeated
  # presentations of the same stimulus, this is y=c(rep(0,times=50),rep(1,times=50))
  # wins - vector indicating NUMBER of time windows (not number of miliseconds 
  # post stimulus onset). For power and phase data, and voltage for comparison
  # with those data, wins <- c(1:166) 
  # ho - number of items in each category in the holdout set (e.g. if ho = 5, 
  # the holdout set includes 5 living and 5 nonliving items)
  # wsz - window size (relevant only for voltage decoding)
  # a - alpha, the elastic net parameter (1 = LASSO, 0 = Ridge)
  # ntpts - number of timepoints present in frequency data (whether you are going
  # to use them all or not!)
  # nf - number of frequencies
  
  # find number of items
	nitems <- dim(d)[1]
	# find number of repetitions of each item. For most analyses, averaging over
	# repeated presentations has already taken place in wavelet.m so this line
	# will set nreps to 1.
	nreps <- nitems/100
	# find number of cross-validation folds. 
	nfolds <- nitems/(ho * 2 * nreps)
	# if voltage data, convert wins from an integer to a vector of window centrepoints
	if (tp=="v"){
	  gap <- 1650/(ntpts-1)
	  wins <- (wins-1)*gap
	}
	
	# initialise outputs
	acc <- matrix(0,nfolds,length(wins))
	ltg <- matrix(0,length(wins),length(wins))
	wwmodels <- list()
	if(return.models) models <- list()

	# Scramble order of living and nonliving. This makes a vector of
	# the numbers 1-100, but 1-50 stay in the first half and 51-100 
	# stay in the second half
	ford <- c(c(1:50)[order(runif(50))], c(51:100)[order(runif(50))])
  # loop through folds
	for(i in c(1:nfolds)){
	  # setup models 
	  ms <- list()
	  # setup fold
    # make a vector of "true"s
	  s <- rep(TRUE, times = nitems/nreps)
		# The first 5 living and the first 5 nonliving stimuli in the scrambled order
	  # are the first holdout set. (The second 5 living and 5 nonliving are the 
	  # next holdout set, etc.)
		hwin <- c(((i-1)*ho + 1):(i*ho))
		hwin <- c(hwin, hwin+50)
		# find the indices of these stimuli within s and replace them with
		# "false"
		hoind <- ford[hwin]
		s[hoind] <- FALSE
		# if there is more than one presentation of each stimulus, hold out all
		# presentations of that stimulus
		s <- rep(s, times = nreps)
		
		for (traint in 1:length(wins)){
		  # get training window 
		  traindata <- get.win(d=d, tp=tp, w = wins[traint], wsz = 50, ntpts = ntpts, nf=nf)
		  # the stimuli marked as "true" go into the training set
		  xtrn <- traindata[s,]
		  ytrn  <- y[s]

		  # train the model on the training set. The measure used is deviance because 
		  # the number of observations per fold is small.
		  m <- cv.glmnet(x = xtrn, y = ytrn, family = "binomial", alpha = a, type.measure = "deviance", lambda = exp(seq(log(0.2), log(0.002), length.out = 100)), parallel = TRUE)
		  if(return.models) ms[[traint]] <- m
		  
		  for (testt in 1:length(wins)){
		    # get data window
		    testdata <- get.win(d=d, tp=tp, w = wins[testt], wsz = 50, ntpts = ntpts, nf=nf)
		    # the stimuli marked as "false" go into the test set
		    xtst <- testdata[!s,]
		    ytst <- y[!s]
		    # give the model a go on the test set
		    p <- assess.glmnet(m, newx = xtst, newy = ytst, s = "lambda.min")
		    # add values to the local temporal generalisation matrix. (This will 
		    # ultimately be divided by the number of folds.)
		    ltg[traint,testt] <- ltg[traint,testt] + p$auc
		    # if the model is being tested on its training window, also fill in the 
		    # matrix of ordinary accuracy
		    if(traint==testt){
		      acc[i,traint] <- p$auc
		      # and print the AUC
		      message(paste("Fold", i, "window", traint, "hold-out AUC:", acc[i,traint]))
		      flush.console()
		    } 
		  }
		# fit one model to the whole of the training window. This will not be used
		# for testing - just for extracting coefficients. We only need to do this
		# once - not every fold - because no data are held out.
		if (i==1) wwmodels[[traint]] <- cv.glmnet(x = traindata, y = y, family="binomial", alpha = 1, type.measure = "deviance")
		# save models
		if(return.models) models[[i]] <- ms
		}
	}
	# divide local temporal generalisation matrix by number of folds
	ltg <- ltg/10
	# output
	if(return.models) o <- list(acc,ltg,wwmodels,models) else o <- list(acc,ltg,wwmodels)
	o
}
		  



