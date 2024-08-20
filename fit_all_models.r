function(path,y=c(rep(0,times=50),rep(1,times=50)),wins=c(1:166),return.models=T){
  # Runs all decoding (all-frequency power, band-specific power, all-frequency
  # phase, band-specific phase, and voltage) and saves output as .Rdata.
  # path - path to /derivatives folder containing .csvs of power, decibel-normalised power and phase for one participant (e.g. "/group/mlr-lab/Saskia/ECoG_LASSO/derivatives/wavelet/sub-01/")
  # y - labels vector. For living/nonliving data, averaged across repeated
  # presentations of the same stimulus, this is y=c(rep(0,times=50),rep(1,times=50))
  # wins - times at which to decode: wins <- seq(1:166)
  # return.models - whether model parameters should be returned
  
  setwd(path)
  # define output path
  outpath <- gsub("/wavelet/","/modelfit/",path)
  dir.create(outpath, recursive=TRUE)
  
  # 1. POWER
  
  # load
  message("Loading power data...")
  d <- read.csv("power.csv",header=F)
  d <- as.matrix(d)
  
  # decode all frequencies
  power <- fit.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:60), return.models=return.models)
  save(power,file=paste(outpath,"power_results.Rdata",sep=""))
  rm(power)
  
  # decode theta (4-7 Hz)
  thetapower <- fit.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:11), return.models=return.models)
  save(thetapower,file=paste(outpath,"theta_power_results.Rdata",sep=""))
  rm(thetapower)
  
  # decode alpha (8-12 Hz)
  alphapower <- fit.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(12:18), return.models=return.models)
  save(alphapower,file=paste(outpath,"alpha_power_results.Rdata",sep=""))
  rm(alphapower)
  
  # decode beta (13-30 Hz)
  betapower <- fit.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(19:31), return.models=return.models)
  save(betapower,file=paste(outpath,"beta_power_results.Rdata",sep=""))
  rm(betapower)
  
  # decode gamma (30-60 Hz)
  gammapower <- fit.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(32:41), return.models=return.models)
  save(gammapower,file=paste(outpath,"gamma_power_results.Rdata",sep=""))
  rm(gammapower)
  
  # decode high gamma (60-200 Hz)
  highgammapower <- fit.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(42:60), return.models=return.models)
  save(highgammapower,file=paste(outpath,"high_gamma_power_results.Rdata",sep=""))
  rm(highgammapower)
  
  rm(d)
  
  # 2. DECIBEL-NORMALISED POWER
  
  # load
  message("Loading decibel-normalised power data...")
  d <- read.csv("dBpower.csv",header=F)
  d <- as.matrix(d)
  
  # decode all frequencies
  dBpower <- fit.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:60), return.models=return.models)
  save(dBpower,file=paste(outpath,"dB_power_results.Rdata",sep=""))
  rm(dBpower)
  
  # decode theta (4-7 Hz)
  dBthetapower <- fit.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:11), return.models=return.models)
  save(dBthetapower,file=paste(outpath,"dB_theta_power_results.Rdata",sep=""))
  rm(dBthetapower)
  
  # decode alpha (8-12 Hz)
  dBalphapower <- fit.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(12:18), return.models=return.models)
  save(dBalphapower,file=paste(outpath,"dB_alpha_power_results.Rdata",sep=""))
  rm(dBalphapower)
  
  # decode beta (13-30 Hz)
  dBbetapower <- fit.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(19:31), return.models=return.models)
  save(dBbetapower,file=paste(outpath,"dB_beta_power_results.Rdata",sep=""))
  rm(dBbetapower)
  
  # decode gamma (30-60 Hz)
  dBgammapower <- fit.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(32:41), return.models=return.models)
  save(dBgammapower,file=paste(outpath,"dB_gamma_power_results.Rdata",sep=""))
  rm(dBgammapower)
  
  # decode high gamma (60-200 Hz)
  dBhighgammapower <- fit.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(42:60), return.models=return.models)
  save(dBhighgammapower,file=paste(outpath,"dB_high_gamma_power_results.Rdata",sep=""))
  rm(dBhighgammapower)
  
  rm(d)
  
  # 3. PHASE
  
  # load
  message("Loading phase data...")
  d <- read.csv("phase.csv",header=F)
  d <- as.matrix(d)
  
  # decode all frequencies
  phase <- fit.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:60), return.models=return.models)
  save(phase,file=paste(outpath,"phase_results.Rdata",sep=""))
  rm(phase)
  
  # decode theta (4-7 Hz)
  thetaphase <- fit.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:11), return.models=return.models)
  save(thetaphase,file=paste(outpath,"theta_phase_results.Rdata",sep=""))
  rm(thetaphase)
  
  # decode alpha (8-12 Hz)
  alphaphase <- fit.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(12:18), return.models=return.models)
  save(alphaphase,file=paste(outpath,"alpha_phase_results.Rdata",sep=""))
  rm(alphaphase)
  
  # decode beta (13-30 Hz)
  betaphase <- fit.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(19:31), return.models=return.models)
  save(betaphase,file=paste(outpath,"beta_phase_results.Rdata",sep=""))
  rm(betaphase)
  
  # decode gamma (30-60 Hz)
  gammaphase <- fit.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(32:41), return.models=return.models)
  save(gammaphase,file=paste(outpath,"gamma_phase_results.Rdata",sep=""))
  rm(gammaphase)
  
  # decode high gamma (60-200 Hz)
  highgammaphase <- fit.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(42:60), return.models=return.models)
  save(highgammaphase,file=paste(outpath,"high_gamma_phase_results.Rdata",sep=""))
  rm(highgammaphase)
  
  rm(d)
  
  # 4. VOLTAGE
  
  # load
  message("Loading voltage data...")
  d <- read.csv("voltage.csv",header=F)
  d <- as.matrix(d)
  
  # decode voltage
  voltage <- fit.models(d, tp="v", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:60), return.models=return.models)
  save(voltage,file=paste(outpath,"voltage_results.Rdata",sep=""))
  rm(voltage)
  
  rm(d)

  message("Done!")
}
