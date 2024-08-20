function(path,y=c(rep(0,times=50),rep(1,times=50)),wins=c(1:166)){
  # Runs all EOG decoding (all-frequency power, band-specific power, all-frequency
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
  message("Loading EOG 'power' data...")
  d <- read.csv("eyes_power.csv",header=F)
  d <- as.matrix(d)
  
  # decode all frequencies
  powereyes <- fit.eog.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:60))
  save(powereyes,file=paste(outpath,"eyes_power_results.Rdata",sep=""))
  rm(powereyes)
  
  # decode theta (4-7 Hz)
  thetapowereyes <- fit.eog.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:11))
  save(thetapowereyes,file=paste(outpath,"eyes_theta_power_results.Rdata",sep=""))
  rm(thetapowereyes)
  
  # decode alpha (8-12 Hz)
  alphapowereyes <- fit.eog.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(12:18))
  save(alphapowereyes,file=paste(outpath,"eyes_alpha_power_results.Rdata",sep=""))
  rm(alphapowereyes)
  
  # decode beta (13-30 Hz)
  betapowereyes <- fit.eog.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(19:31))
  save(betapowereyes,file=paste(outpath,"eyes_beta_power_results.Rdata",sep=""))
  rm(betapowereyes)
  
  # decode gamma (30-60 Hz)
  gammapowereyes <- fit.eog.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(32:41))
  save(gammapowereyes,file=paste(outpath,"eyes_gamma_power_results.Rdata",sep=""))
  rm(gammapowereyes)
  
  # decode high gamma (60-200 Hz)
  highgammapowereyes <- fit.eog.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(42:60))
  save(highgammapowereyes,file=paste(outpath,"eyes_high_gamma_power_results.Rdata",sep=""))
  rm(highgammapowereyes)
  
  rm(d)
  
  # 2. DECIBEL-NORMALISED POWER
  
  # load
  message("Loading EOG 'decibel-normalised power' data...")
  d <- read.csv("eyes_dBpower.csv",header=F)
  d <- as.matrix(d)
  
  # decode all frequencies
  dBpowereyes <- fit.eog.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:60))
  save(dBpowereyes,file=paste(outpath,"eyes_dB_power_results.Rdata",sep=""))
  rm(dBpowereyes)
  
  # decode theta (4-7 Hz)
  dBthetapowereyes <- fit.eog.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:11))
  save(dBthetapowereyes,file=paste(outpath,"eyes_dB_theta_power_results.Rdata",sep=""))
  rm(dBthetapowereyes)
  
  # decode alpha (8-12 Hz)
  dBalphapowereyes <- fit.eog.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(12:18))
  save(dBalphapowereyes,file=paste(outpath,"eyes_dB_alpha_power_results.Rdata",sep=""))
  rm(dBalphapowereyes)
  
  # decode beta (13-30 Hz)
  dBbetapowereyes <- fit.eog.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(19:31))
  save(dBbetapowereyes,file=paste(outpath,"eyes_dB_beta_power_results.Rdata",sep=""))
  rm(dBbetapowereyes)
  
  # decode gamma (30-60 Hz)
  dBgammapowereyes <- fit.eog.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(32:41))
  save(dBgammapowereyes,file=paste(outpath,"eyes_dB_gamma_power_results.Rdata",sep=""))
  rm(dBgammapowereyes)
  
  # decode high gamma (60-200 Hz)
  dBhighgammapowereyes <- fit.eog.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(42:60))
  save(dBhighgammapowereyes,file=paste(outpath,"eyes_dB_high_gamma_power_results.Rdata",sep=""))
  rm(dBhighgammapowereyes)
  
  rm(d)
  
  # 3. PHASE
 
  # load
  message("Loading EOG 'phase' data...")
  d <- read.csv("eyes_phase.csv",header=F)
  d <- as.matrix(d)
  
  # decode all frequencies
  phaseeyes <- fit.eog.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:60))
  save(phaseeyes,file=paste(outpath,"eyes_phase_results.Rdata",sep=""))
  rm(phaseeyes)
  
  # decode theta (4-7 Hz)
  thetaphaseeyes <- fit.eog.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:11))
  save(thetaphaseeyes,file=paste(outpath,"eyes_theta_phase_results.Rdata",sep=""))
  rm(thetaphaseeyes)
  
  # decode alpha (8-12 Hz)
  alphaphaseeyes <- fit.eog.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(12:18))
  save(alphaphaseeyes,file=paste(outpath,"eyes_alpha_phase_results.Rdata",sep=""))
  rm(alphaphaseeyes)
  
  # decode beta (13-30 Hz)
  betaphaseeyes <- fit.eog.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(19:31))
  save(betaphaseeyes,file=paste(outpath,"eyes_beta_phase_results.Rdata",sep=""))
  rm(betaphaseeyes)
  
  # decode gamma (30-60 Hz)
  gammaphaseeyes <- fit.eog.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(32:41))
  save(gammaphaseeyes,file=paste(outpath,"eyes_gamma_phase_results.Rdata",sep=""))
  rm(gammaphaseeyes)
  
  # decode high gamma (60-200 Hz)
  highgammaphaseeyes <- fit.eog.models(d, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(42:60))
  save(highgammaphaseeyes,file=paste(outpath,"eyes_high_gamma_phase_results.Rdata",sep=""))
  rm(highgammaphaseeyes)
  
  rm(d)
  
  # 4. VOLTAGE
  
  # load
  message("Loading EOG 'voltage' data...")
  d <- read.csv("eyes_voltage.csv",header=F)
  d <- as.matrix(d)
  
  # decode voltage
  voltageeyes <- fit.eog.models(d, tp="v", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:60))
  save(voltageeyes,file=paste(outpath,"eyes_voltage_results.Rdata",sep=""))
  rm(voltageeyes)
  
  rm(d)

  message("Done!")
}
