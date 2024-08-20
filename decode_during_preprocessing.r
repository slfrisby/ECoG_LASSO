function(path,y=c(rep(0,times=50),rep(1,times=50)),wins=c(1:166)){
  # Runs decoding during preprocessing (all-frequency power, band-specific power, all-frequency
  # phase, band-specific phase, and voltage) and saves output as .Rdata.
  # path - path to /derivatives folder containing .csvs of power, decibel-normalised power and phase for one participant (e.g. "/group/mlr-lab/Saskia/ECoG_LASSO/derivatives/wavelet/sub-01/")
  # y - labels vector. For living/nonliving data, averaged across repeated
  # presentations of the same stimulus, this is y=c(rep(0,times=50),rep(1,times=50))
  # wins - times at which to decode: wins <- seq(1:166)
  
  setwd(path)
  # define output path
  outpath <- gsub("/wavelet/","/modelfit/",path)
  dir.create(outpath, recursive=TRUE)
  
  # # 1. POWER
  # 
  # # load
  # message("Loading power data at each preprocessing stage...")
  # raw <- read.csv("raw_power.csv",header=F)
  # raw <- as.matrix(raw)
  # filtered <- read.csv("filtered_power.csv",header=F)
  # filtered <- as.matrix(filtered)
  # channelsrejected <- read.csv("channelsrejected_power.csv",header=F)
  # channelsrejected <- as.matrix(channelsrejected)
  # rereferenced <- read.csv("rereferenced_power.csv",header=F)
  # rereferenced <- as.matrix(rereferenced)
  # trialsrejected <- read.csv("trialsrejected_power.csv",header=F)
  # trialsrejected <- as.matrix(trialsrejected)
  # 
  # # decode all frequencies
  # powerddp <- fit.models.during.preprocessing(raw=raw, filtered=filtered, channelsrejected=channelsrejected, rereferenced=rereferenced, trialsrejected=trialsrejected, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:60))
  # save(powerddp,file=paste(outpath,"ddp_power_results.Rdata",sep=""))
  # rm(powerddp)
  # 
  # # decode theta (4-7 Hz)
  # thetapowerddp <- fit.models.during.preprocessing(raw=raw, filtered=filtered, channelsrejected=channelsrejected, rereferenced=rereferenced, trialsrejected=trialsrejected, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:11))
  # save(thetapowerddp,file=paste(outpath,"ddp_theta_power_results.Rdata",sep=""))
  # rm(thetapowerddp)
  # 
  # # decode alpha (8-12 Hz)
  # alphapowerddp <- fit.models.during.preprocessing(raw=raw, filtered=filtered, channelsrejected=channelsrejected, rereferenced=rereferenced, trialsrejected=trialsrejected, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(12:18))
  # save(alphapowerddp,file=paste(outpath,"ddp_alpha_power_results.Rdata",sep=""))
  # rm(alphapowerddp)
  # 
  # # decode beta (13-30 Hz)
  # betapowerddp <- fit.models.during.preprocessing(raw=raw, filtered=filtered, channelsrejected=channelsrejected, rereferenced=rereferenced, trialsrejected=trialsrejected, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(19:31))
  # save(betapowerddp,file=paste(outpath,"ddp_beta_power_results.Rdata",sep=""))
  # rm(betapowerddp)
  # 
  # # decode gamma (30-60 Hz)
  # gammapowerddp <- fit.models.during.preprocessing(raw=raw, filtered=filtered, channelsrejected=channelsrejected, rereferenced=rereferenced, trialsrejected=trialsrejected, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(32:41))
  # save(gammapowerddp,file=paste(outpath,"ddp_gamma_power_results.Rdata",sep=""))
  # rm(gammapowerddp)
  # 
  # # decode high gamma (60-200 Hz)
  # highgammapowerddp <- fit.models.during.preprocessing(raw=raw, filtered=filtered, channelsrejected=channelsrejected, rereferenced=rereferenced, trialsrejected=trialsrejected, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(42:60))
  # save(highgammapowerddp,file=paste(outpath,"ddp_high_gamma_power_results.Rdata",sep=""))
  # rm(highgammapowerddp)
  # 
  # rm(raw,filtered,channelsrejected,rereferenced,trialsrejected)
  
  # 2. DECIBEL-NORMALISED POWER
  
  # load
  message("Loading decibel-normalised power data at each preprocessing stage...")
  raw <- read.csv("raw_dBpower.csv",header=F)
  raw <- as.matrix(raw)
  filtered <- read.csv("filtered_dBpower.csv",header=F)
  filtered <- as.matrix(filtered)
  channelsrejected <- read.csv("channelsrejected_dBpower.csv",header=F)
  channelsrejected <- as.matrix(channelsrejected)
  rereferenced <- read.csv("rereferenced_dBpower.csv",header=F)
  rereferenced <- as.matrix(rereferenced)
  trialsrejected <- read.csv("trialsrejected_dBpower.csv",header=F)
  trialsrejected <- as.matrix(trialsrejected)
  
  # decode all frequencies
  dBpowerddp <- fit.models.during.preprocessing(raw=raw, filtered=filtered, channelsrejected=channelsrejected, rereferenced=rereferenced, trialsrejected=trialsrejected, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:60))
  save(dBpowerddp,file=paste(outpath,"ddp_dB_power_results.Rdata",sep=""))
  rm(dBpowerddp)
  
  # decode theta (4-7 Hz)
  dBthetapowerddp <- fit.models.during.preprocessing(raw=raw, filtered=filtered, channelsrejected=channelsrejected, rereferenced=rereferenced, trialsrejected=trialsrejected, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:11))
  save(dBthetapowerddp,file=paste(outpath,"ddp_dB_theta_power_results.Rdata",sep=""))
  rm(dBthetapowerddp)
  
  # decode alpha (8-12 Hz)
  dBalphapowerddp <- fit.models.during.preprocessing(raw=raw, filtered=filtered, channelsrejected=channelsrejected, rereferenced=rereferenced, trialsrejected=trialsrejected, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(12:18))
  save(dBalphapowerddp,file=paste(outpath,"ddp_dB_alpha_power_results.Rdata",sep=""))
  rm(dBalphapowerddp)
  
  # decode beta (13-30 Hz)
  dBbetapowerddp <- fit.models.during.preprocessing(raw=raw, filtered=filtered, channelsrejected=channelsrejected, rereferenced=rereferenced, trialsrejected=trialsrejected, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(19:31))
  save(dBbetapowerddp,file=paste(outpath,"ddp_dB_beta_power_results.Rdata",sep=""))
  rm(dBbetapowerddp)
  
  # decode gamma (30-60 Hz)
  dBgammapowerddp <- fit.models.during.preprocessing(raw=raw, filtered=filtered, channelsrejected=channelsrejected, rereferenced=rereferenced, trialsrejected=trialsrejected, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(32:41))
  save(dBgammapowerddp,file=paste(outpath,"ddp_dB_gamma_power_results.Rdata",sep=""))
  rm(dBgammapowerddp)
  
  # decode high gamma (60-200 Hz)
  dBhighgammapowerddp <- fit.models.during.preprocessing(raw=raw, filtered=filtered, channelsrejected=channelsrejected, rereferenced=rereferenced, trialsrejected=trialsrejected, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(42:60))
  save(dBhighgammapowerddp,file=paste(outpath,"ddp_dB_high_gamma_power_results.Rdata",sep=""))
  rm(dBhighgammapowerddp)
  
  rm(raw,filtered,channelsrejected,rereferenced,trialsrejected)
  
  # 3. PHASE
  
  # load
  message("Loading phase data at each preprocessing stage...")
  raw <- read.csv("raw_phase.csv",header=F)
  raw <- as.matrix(raw)
  filtered <- read.csv("filtered_phase.csv",header=F)
  filtered <- as.matrix(filtered)
  channelsrejected <- read.csv("channelsrejected_phase.csv",header=F)
  channelsrejected <- as.matrix(channelsrejected)
  rereferenced <- read.csv("rereferenced_phase.csv",header=F)
  rereferenced <- as.matrix(rereferenced)
  trialsrejected <- read.csv("trialsrejected_phase.csv",header=F)
  trialsrejected <- as.matrix(trialsrejected)
  
  # decode all frequencies
  phaseddp <- fit.models.during.preprocessing(raw=raw, filtered=filtered, channelsrejected=channelsrejected, rereferenced=rereferenced, trialsrejected=trialsrejected, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:60))
  save(phaseddp,file=paste(outpath,"ddp_phase_results.Rdata",sep=""))
  rm(phaseddp)
  
  # decode theta (4-7 Hz)
  thetaphaseddp <- fit.models.during.preprocessing(raw=raw, filtered=filtered, channelsrejected=channelsrejected, rereferenced=rereferenced, trialsrejected=trialsrejected, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:11))
  save(thetaphaseddp,file=paste(outpath,"ddp_theta_phase_results.Rdata",sep=""))
  rm(thetaphaseddp)
  
  # decode alpha (8-12 Hz)
  alphaphaseddp <- fit.models.during.preprocessing(raw=raw, filtered=filtered, channelsrejected=channelsrejected, rereferenced=rereferenced, trialsrejected=trialsrejected, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(12:18))
  save(alphaphaseddp,file=paste(outpath,"ddp_alpha_phase_results.Rdata",sep=""))
  rm(alphaphaseddp)
  
  # decode beta (13-30 Hz)
  betaphaseddp <- fit.models.during.preprocessing(raw=raw, filtered=filtered, channelsrejected=channelsrejected, rereferenced=rereferenced, trialsrejected=trialsrejected, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(19:31))
  save(betaphaseddp,file=paste(outpath,"ddp_beta_phase_results.Rdata",sep=""))
  rm(betaphaseddp)
  
  # decode gamma (30-60 Hz)
  gammaphaseddp <- fit.models.during.preprocessing(raw=raw, filtered=filtered, channelsrejected=channelsrejected, rereferenced=rereferenced, trialsrejected=trialsrejected, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(32:41))
  save(gammaphaseddp,file=paste(outpath,"ddp_gamma_phase_results.Rdata",sep=""))
  rm(gammaphaseddp)
  
  # decode high gamma (60-200 Hz)
  highgammaphaseddp <- fit.models.during.preprocessing(raw=raw, filtered=filtered, channelsrejected=channelsrejected, rereferenced=rereferenced, trialsrejected=trialsrejected, tp="f", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(42:60))
  save(highgammaphaseddp,file=paste(outpath,"ddp_high_gamma_phase_results.Rdata",sep=""))
  rm(highgammaphaseddp)
  
  rm(raw,filtered,channelsrejected,rereferenced,trialsrejected)
  
  # 4. VOLTAGE
  
  # load
  message("Loading voltage data at each preprocessing stage...")
  raw <- read.csv("raw_voltage.csv",header=F)
  raw <- as.matrix(raw)
  filtered <- read.csv("filtered_voltage.csv",header=F)
  filtered <- as.matrix(filtered)
  channelsrejected <- read.csv("channelsrejected_voltage.csv",header=F)
  channelsrejected <- as.matrix(channelsrejected)
  rereferenced <- read.csv("rereferenced_voltage.csv",header=F)
  rereferenced <- as.matrix(rereferenced)
  trialsrejected <- read.csv("trialsrejected_voltage.csv",header=F)
  trialsrejected <- as.matrix(trialsrejected)
  
  # decode voltage
  voltageddp <- fit.models.during.preprocessing(raw=raw, filtered=filtered, channelsrejected=channelsrejected, rereferenced=rereferenced, trialsrejected=trialsrejected, tp="v", y=y, wins=wins, ho = 5, wsz = 50, a = 1, ntpts = 166, nf = c(1:60))
  save(voltageddp,file=paste(outpath,"ddp_voltage_results.Rdata",sep=""))
  rm(voltageddp)
  
  rm(raw,filtered,channelsrejected,rereferenced,trialsrejected)
  
  message("Done!")
}
