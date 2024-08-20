function(path,ntpts=166,wwidth=50,fileformat="png"){
  # Takes results of fit.all.models, plots figures for individual participants,
  # collates results across participants, and plots figures of group results.
  # path - path to /modelfit directory of results (path <- "Z:/Saskia/ECoG_LASSO/derivatives/modelfit/")
  # ntpts - number of time points
  # wwidth - width of window in ms
  # fileformat - format to save pictures. Options: "png" (default),"svg"
  
  # get dependencies
  library(glmnet)
  library(fields)
  library(effsize)
  library(lsa)
  library(segmented)
  library(flux)
  setwd("Z:/Saskia/ECoG_LASSO/scripts")
  plot.cis <- dget("plot_cis.r")
  stars <- dget("stars.r")
  plot.coefs <- dget("plot_coefs.r")
  get.clust <- dget("get_clust.r")
  plot.clacc <- dget("plot_clacc.r")
  plot.cluster.compare <- dget("plot_cluster_compare.r")
  plot.auc.widen <- dget("plot_auc_widen.r")
  setwd(path)

  # make output directory for collated results
  grouppath <- paste(path,"/group/",sep="")
  dir.create(grouppath, recursive=TRUE)
  # make output directories for figures
  figurespath <- gsub("modelfit","figures",path)
  groupfigurespath <- paste(figurespath,"/group",sep="")
  dir.create(groupfigurespath, recursive=TRUE)
  # make work directory for coefficients
  coefwork <- gsub("derivatives/modelfit","work/coefficients/R",path)
  dir.create(coefwork, recursive=TRUE)
  
  # set participant numbers, excluding bad participants
  subs <- c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","20","21","22")
  # only some participants have EOG
  eogsubs <- c("01","02","03","06","22")
  
  # initialise a lot of variables!
  dBgrouppowerrawacc <- matrix(0,length(subs),ntpts)
  dBgrouppowerfilteredacc <- matrix(0,length(subs),ntpts)
  dBgrouppowerchannelsrejectedacc <- matrix(0,length(subs),ntpts)
  dBgrouppowerrereferencedacc <- matrix(0,length(subs),ntpts)
  dBgrouppowertrialsrejectedacc <- matrix(0,length(subs),ntpts)
  groupphaserawacc <- matrix(0,length(subs),ntpts)
  groupphasefilteredacc <- matrix(0,length(subs),ntpts)
  groupphasechannelsrejectedacc <- matrix(0,length(subs),ntpts)
  groupphaserereferencedacc <- matrix(0,length(subs),ntpts)
  groupphasetrialsrejectedacc <- matrix(0,length(subs),ntpts)
  groupvoltagerawacc <- matrix(0,length(subs),ntpts)
  groupvoltagefilteredacc <- matrix(0,length(subs),ntpts)
  groupvoltagechannelsrejectedacc <- matrix(0,length(subs),ntpts)
  groupvoltagerereferencedacc <- matrix(0,length(subs),ntpts)
  groupvoltagetrialsrejectedacc <- matrix(0,length(subs),ntpts)
  dBgrouppoweracc <- matrix(0,length(subs),ntpts)
  dBgroupthetapoweracc <- matrix(0,length(subs),ntpts)
  dBgroupalphapoweracc <- matrix(0,length(subs),ntpts)
  dBgroupbetapoweracc <- matrix(0,length(subs),ntpts)
  dBgroupgammapoweracc <- matrix(0,length(subs),ntpts)
  dBgrouphighgammapoweracc <- matrix(0,length(subs),ntpts)
  dBallpowerltg <- array(0,dim=c(ntpts,ntpts,length(subs)))
  dBallthetapowerltg <- array(0,dim=c(ntpts,ntpts,length(subs)))
  dBallalphapowerltg <- array(0,dim=c(ntpts,ntpts,length(subs)))
  dBallbetapowerltg <- array(0,dim=c(ntpts,ntpts,length(subs)))
  dBallgammapowerltg <- array(0,dim=c(ntpts,ntpts,length(subs)))
  dBallhighgammapowerltg <- array(0,dim=c(ntpts,ntpts,length(subs)))
  groupphaseacc <- matrix(0,length(subs),ntpts)
  groupthetaphaseacc <- matrix(0,length(subs),ntpts)
  groupalphaphaseacc <- matrix(0,length(subs),ntpts)
  groupbetaphaseacc <- matrix(0,length(subs),ntpts)
  groupgammaphaseacc <- matrix(0,length(subs),ntpts)
  grouphighgammaphaseacc <- matrix(0,length(subs),ntpts)
  allphaseltg <- array(0,dim=c(ntpts,ntpts,length(subs)))
  allthetaphaseltg <- array(0,dim=c(ntpts,ntpts,length(subs)))
  allalphaphaseltg <- array(0,dim=c(ntpts,ntpts,length(subs)))
  allbetaphaseltg <- array(0,dim=c(ntpts,ntpts,length(subs)))
  allgammaphaseltg <- array(0,dim=c(ntpts,ntpts,length(subs)))
  allhighgammaphaseltg <- array(0,dim=c(ntpts,ntpts,length(subs)))
  groupvoltageacc <- matrix(0,length(subs),ntpts)
  allvoltageltg <- array(0,dim=c(ntpts,ntpts,length(subs)))
  # only 5 participants have EOG data
  dBgrouppowereyes <- matrix(0,5,ntpts)
  dBgroupthetapowereyes <- matrix(0,5,ntpts)
  dBgroupalphapowereyes <- matrix(0,5,ntpts)
  dBgroupbetapowereyes <- matrix(0,5,ntpts)
  dBgroupgammapowereyes <- matrix(0,5,ntpts)
  dBgrouphighgammapowereyes <- matrix(0,5,ntpts)
  groupphaseeyes <- matrix(0,5,ntpts)
  groupthetaphaseeyes <- matrix(0,5,ntpts)
  groupalphaphaseeyes <- matrix(0,5,ntpts)
  groupbetaphaseeyes <- matrix(0,5,ntpts)
  groupgammaphaseeyes <- matrix(0,5,ntpts)
  grouphighgammaphaseeyes <- matrix(0,5,ntpts)
  groupvoltageeyes <- matrix(0,5,ntpts)

  for(i in c(1:length(subs))){
    
    # make output directory for figures
    fpath <- paste(figurespath,"/sub-",subs[i],sep="")
    dir.create(fpath, recursive=TRUE)
    # get work directory
    work <- gsub("derivatives/modelfit",paste("work/sub-",subs[i],sep=""),path)

    # 1. PREPROCESSING CHECKS

    # decibel-normalised power

    load(paste(path,"/sub-",subs[i],"/ddp_dB_power_results.Rdata",sep=""))
    rawacc <- dBpowerddp[[1]]
    filteredacc <- dBpowerddp[[2]]
    channelsrejectedacc <- dBpowerddp[[3]]
    rereferencedacc <- dBpowerddp[[4]]
    trialsrejectedacc <- dBpowerddp[[5]]
    rm(dBpowerddp)
    if (fileformat=="png"){png(paste(fpath,"/dB_DDP_power.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(fpath,"/dB_DDP_power.svg",sep=""))}
    plot.cis(rawacc, colour="#583C88")
    stars(d1=rawacc,colour1="#583C88",yval=1)
    plot.cis(filteredacc, colour="#9B4FA5", newflag = F)
    stars(d1=filteredacc,colour1="#9B4FA5",yval=0.975)
    plot.cis(channelsrejectedacc, colour="#C4619E", newflag = F)
    stars(d1=channelsrejectedacc,colour1="#C4619E",yval=0.95)
    plot.cis(rereferencedacc, colour="#DE89B0", newflag = F)
    stars(d1=rereferencedacc,colour1="#DE89B0", yval=0.925)
    plot.cis(trialsrejectedacc, colour="#F7BCCB", newflag = F)
    stars(d1=trialsrejectedacc,colour1="#F7BCCB",yval=0.9)
    title(paste("Sub-",subs[i],sep=""),cex.main = 2)
    dev.off()
    # collate results
    dBgrouppowerrawacc[i,] <- colMeans(rawacc)
    dBgrouppowerfilteredacc[i,] <- colMeans(filteredacc)
    dBgrouppowerchannelsrejectedacc[i,] <- colMeans(channelsrejectedacc)
    dBgrouppowerrereferencedacc[i,] <- colMeans(rereferencedacc)
    dBgrouppowertrialsrejectedacc[i,] <- colMeans(trialsrejectedacc)
    rm(rawacc,filteredacc,channelsrejectedacc,rereferencedacc,trialsrejectedacc)
    
    # phase

    load(paste(path,"/sub-",subs[i],"/ddp_phase_results.Rdata",sep=""))
    rawacc <- phaseddp[[1]]
    filteredacc <- phaseddp[[2]]
    channelsrejectedacc <- phaseddp[[3]]
    rereferencedacc <- phaseddp[[4]]
    trialsrejectedacc <- phaseddp[[5]]
    rm(phaseddp)
    if (fileformat=="png"){png(paste(fpath,"/DDP_phase.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(fpath,"/DDP_phase.svg",sep=""))}
    plot.cis(rawacc, colour="#583C88")
    stars(d1=rawacc,colour1="#583C88",yval=1)
    plot.cis(filteredacc, colour="#9B4FA5", newflag = F)
    stars(d1=filteredacc,colour1="#9B4FA5",yval=0.975)
    plot.cis(channelsrejectedacc, colour="#C4619E", newflag = F)
    stars(d1=channelsrejectedacc,colour1="#C4619E",yval=0.95)
    plot.cis(rereferencedacc, colour="#DE89B0", newflag = F)
    stars(d1=rereferencedacc,colour1="#DE89B0", yval=0.925)
    plot.cis(trialsrejectedacc, colour="#F7BCCB", newflag = F)
    stars(d1=trialsrejectedacc,colour1="#F7BCCB",yval=0.9)
    title(paste("Sub-",subs[i],sep=""),cex.main = 2)
    dev.off()
    # collate results
    groupphaserawacc[i,] <- colMeans(rawacc)
    groupphasefilteredacc[i,] <- colMeans(filteredacc)
    groupphasechannelsrejectedacc[i,] <- colMeans(channelsrejectedacc)
    groupphaserereferencedacc[i,] <- colMeans(rereferencedacc)
    groupphasetrialsrejectedacc[i,] <- colMeans(trialsrejectedacc)
    rm(rawacc,filteredacc,channelsrejectedacc,rereferencedacc,trialsrejectedacc)

    # voltage

    load(paste(path,"/sub-",subs[i],"/ddp_voltage_results.Rdata",sep=""))
    rawacc <- voltageddp[[1]]
    filteredacc <- voltageddp[[2]]
    channelsrejectedacc <- voltageddp[[3]]
    rereferencedacc <- voltageddp[[4]]
    trialsrejectedacc <- voltageddp[[5]]
    rm(voltageddp)
    if (fileformat=="png"){png(paste(fpath,"/DDP_voltage.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(fpath,"/DDP_voltage.svg",sep=""))}
    plot.cis(rawacc, colour="#583C88")
    stars(d1=rawacc,colour1="#583C88",yval=1)
    plot.cis(filteredacc, colour="#9B4FA5", newflag = F)
    stars(d1=filteredacc,colour1="#9B4FA5",yval=0.975)
    plot.cis(channelsrejectedacc, colour="#C4619E", newflag = F)
    stars(d1=channelsrejectedacc,colour1="#C4619E",yval=0.95)
    plot.cis(rereferencedacc, colour="#DE89B0", newflag = F)
    stars(d1=rereferencedacc,colour1="#DE89B0", yval=0.925)
    plot.cis(trialsrejectedacc, colour="#F7BCCB", newflag = F)
    stars(d1=trialsrejectedacc,colour1="#F7BCCB",yval=0.9)
    title(paste("Sub-",subs[i],sep=""),cex.main = 2)
    dev.off()
    # collate results
    groupvoltagerawacc[i,] <- colMeans(rawacc)
    groupvoltagefilteredacc[i,] <- colMeans(filteredacc)
    groupvoltagechannelsrejectedacc[i,] <- colMeans(channelsrejectedacc)
    groupvoltagerereferencedacc[i,] <- colMeans(rereferencedacc)
    groupvoltagetrialsrejectedacc[i,] <- colMeans(trialsrejectedacc)
    rm(rawacc,filteredacc,channelsrejectedacc,rereferencedacc,trialsrejectedacc)

    # 2. MAIN ANALYSIS 
    
    # load voltage - accuracy, temporal generalisation, and models fit to all of
    # the data in each window
    load(paste(path,"/sub-",subs[i],"/voltage_results.Rdata",sep=""))
    voltageacc <- voltage[[1]]
    voltageltg <- voltage[[2]]
    voltagemodels <- voltage[[3]]
    rm(voltage)
    
    # power
    
    # decibel-normalised power

    # load all decibel-normalised power data
    
    load(paste(path,"/sub-",subs[i],"/dB_power_results.Rdata",sep=""))
    dBpoweracc <- dBpower[[1]]
    dBpowerltg <- dBpower[[2]]
    dBpowermodels <- dBpower[[3]]
    load(paste(path,"/sub-",subs[i],"/dB_theta_power_results.Rdata",sep=""))
    dBthetapoweracc <- dBthetapower[[1]]
    dBthetapowerltg <- dBthetapower[[2]]
    dBthetapowermodels <- dBthetapower[[3]]
    load(paste(path,"/sub-",subs[i],"/dB_alpha_power_results.Rdata",sep=""))
    dBalphapoweracc <- dBalphapower[[1]]
    dBalphapowerltg <- dBalphapower[[2]]
    dBalphapowermodels <- dBalphapower[[3]]
    load(paste(path,"/sub-",subs[i],"/dB_beta_power_results.Rdata",sep=""))
    dBbetapoweracc <- dBbetapower[[1]]
    dBbetapowerltg <- dBbetapower[[2]]
    dBbetapowermodels <- dBbetapower[[3]]
    load(paste(path,"/sub-",subs[i],"/dB_gamma_power_results.Rdata",sep=""))
    dBgammapoweracc <- dBgammapower[[1]]
    dBgammapowerltg <- dBgammapower[[2]]
    dBgammapowermodels <- dBgammapower[[3]]
    load(paste(path,"/sub-",subs[i],"/dB_high_gamma_power_results.Rdata",sep=""))
    dBhighgammapoweracc <- dBhighgammapower[[1]]
    dBhighgammapowerltg <- dBhighgammapower[[2]]
    dBhighgammapowermodels <- dBhighgammapower[[3]]
    rm(dBpower,dBthetapower,dBalphapower,dBbetapower,dBgammapower,dBhighgammapower)
    
    # plot all frequencies against voltage
    if (fileformat=="png"){png(paste(fpath,"/dB_power_freq_vs_volt.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(fpath,"/dB_power_freq_vs_volt.svg",sep=""))}
    plot.cis(voltageacc, colour="#996035")
    plot.cis(dBpoweracc, colour="#0066FF", newflag = F)
    title(paste("Sub-",subs[i],sep=""), cex.main = 2)
    # add stars to indicate difference from chance
    stars(d1=voltageacc,colour1="#996035",yval=1)
    stars(d1=dBpoweracc,colour1="#0066FF",yval=0.975)
    # add stars to indicate difference from each other
    stars(d1=voltageacc,d2=dBpoweracc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
    #save
    dev.off()
    
    # collate results
    dBgrouppoweracc[i,] <- colMeans(dBpoweracc)
    rm(dBpoweracc)
    
    # plot all bands
    if (fileformat=="png"){png(paste(fpath,"/dB_power_bands.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(fpath,"/dB_power_bands.svg",sep=""))}
    plot.cis(dBthetapoweracc, colour="#A50026")
    stars(d1=dBthetapoweracc,colour1="#A50026",yval=1)
    plot.cis(dBalphapoweracc, colour="#F46D43", newflag = F)
    stars(d1=dBalphapoweracc,colour1="#F46D43",yval=0.975)
    plot.cis(dBbetapoweracc, colour="#FDCC3F", newflag = F)
    stars(d1=dBbetapoweracc,colour1="#FDCC3F",yval=0.95)
    plot.cis(dBgammapoweracc, colour="#66BD63", newflag = F)
    stars(d1=dBgammapoweracc,colour1="#66BD63",yval=0.925)
    plot.cis(dBhighgammapoweracc, colour="#006837", newflag = F)
    stars(d1=dBhighgammapoweracc,colour1="#006837",yval=0.9)
    title(paste("Sub-",subs[i],sep=""), cex.main = 2)
    dev.off()
    
    # collate results
    dBgroupthetapoweracc[i,] <- colMeans(dBthetapoweracc)
    dBgroupalphapoweracc[i,] <- colMeans(dBalphapoweracc)
    dBgroupbetapoweracc[i,] <- colMeans(dBbetapoweracc)
    dBgroupgammapoweracc[i,] <- colMeans(dBgammapoweracc)
    dBgrouphighgammapoweracc[i,] <- colMeans(dBhighgammapoweracc)
    rm(dBthetapoweracc,dBalphapoweracc,dBbetapoweracc,dBgammapoweracc,dBhighgammapoweracc)

    
    # plot local temporal generalisation matrices. Reversing the index means that the axes are oriented correctly
    if (fileformat=="png"){png(paste(fpath,"/dB_power_LTG.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(fpath,"/power_LTG",sep=""))}
    image.plot(seq(0,1650,by=10),seq(0,1650,by=10),dBpowerltg[,ntpts:1],zlim=c(0.4,0.8),xaxt="n",yaxt="n",xlab="Test time (ms) →",ylab="← Train time (ms)",cex.lab = 1.5,axis.args=list(cex.axis=1.5))
    box()
    axis(1, at = seq(0,1650,by=500),labels= c(0,500,1000,1500),cex.axis = 1.5)
    axis(2, at = seq(1650,0,by=-500),labels= c(0,500,1000,1500),cex.axis = 1.5)
    title(paste("Sub-",subs[i]," - all frequencies",sep=""),cex.main = 2)
    dev.off()
    if (fileformat=="png"){png(paste(fpath,"/dB_theta_power_LTG.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(fpath,"/dB_theta_power_LTG",sep=""))}
    image.plot(seq(0,1650,by=10),seq(0,1650,by=10),dBthetapowerltg[,ntpts:1],zlim=c(0.4,0.8),xaxt="n",yaxt="n",xlab="Test time (ms) →",ylab="← Train time (ms)",cex.lab = 1.5,axis.args=list(cex.axis=1.5))
    box()
    axis(1, at = seq(0,1650,by=500),labels= c(0,500,1000,1500),cex.axis = 1.5)
    axis(2, at = seq(1650,0,by=-500),labels= c(0,500,1000,1500),cex.axis = 1.5)
    title(paste("Sub-",subs[i]," - theta",sep=""),cex.main = 2)
    dev.off()
    if (fileformat=="png"){png(paste(fpath,"/dB_alpha_power_LTG.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(fpath,"/dB_alpha_power_LTG",sep=""))}
    image.plot(seq(0,1650,by=10),seq(0,1650,by=10),dBalphapowerltg[,ntpts:1],zlim=c(0.4,0.8),xaxt="n",yaxt="n",xlab="Test time (ms) →",ylab="← Train time (ms)",cex.lab = 1.5,axis.args=list(cex.axis=1.5))
    box()
    axis(1, at = seq(0,1650,by=500),labels= c(0,500,1000,1500),cex.axis = 1.5)
    axis(2, at = seq(1650,0,by=-500),labels= c(0,500,1000,1500),cex.axis = 1.5)
    title(paste("Sub-",subs[i]," - alpha",sep=""),cex.main = 2)
    dev.off()
    if (fileformat=="png"){png(paste(fpath,"/dB_beta_power_LTG.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(fpath,"/dB_beta_power_LTG",sep=""))}
    image.plot(seq(0,1650,by=10),seq(0,1650,by=10),dBbetapowerltg[,ntpts:1],zlim=c(0.4,0.8),xaxt="n",yaxt="n",xlab="Test time (ms) →",ylab="← Train time (ms)",cex.lab = 1.5,axis.args=list(cex.axis=1.5))
    box()
    axis(1, at = seq(0,1650,by=500),labels= c(0,500,1000,1500),cex.axis = 1.5)
    axis(2, at = seq(1650,0,by=-500),labels= c(0,500,1000,1500),cex.axis = 1.5)
    title(paste("Sub-",subs[i]," - beta",sep=""),cex.main = 2)
    dev.off()
    if (fileformat=="png"){png(paste(fpath,"/dB_gamma_power_LTG.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(fpath,"/dB_gamma_power_LTG",sep=""))}
    image.plot(seq(0,1650,by=10),seq(0,1650,by=10),dBgammapowerltg[,ntpts:1],zlim=c(0.4,0.8),xaxt="n",yaxt="n",xlab="Test time (ms) →",ylab="← Train time (ms)",cex.lab = 1.5,axis.args=list(cex.axis=1.5))
    box()
    axis(1, at = seq(0,1650,by=500),labels= c(0,500,1000,1500),cex.axis = 1.5)
    axis(2, at = seq(1650,0,by=-500),labels= c(0,500,1000,1500),cex.axis = 1.5)
    title(paste("Sub-",subs[i]," - gamma",sep=""),cex.main = 2)
    dev.off()
    if (fileformat=="png"){png(paste(fpath,"/dB_high_gamma_power_LTG.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(fpath,"/dB_high_gamma_power_LTG",sep=""))}
    image.plot(seq(0,1650,by=10),seq(0,1650,by=10),dBhighgammapowerltg[,ntpts:1],zlim=c(0.4,0.8),xaxt="n",yaxt="n",xlab="Test time (ms) →",ylab="← Train time (ms)",cex.lab = 1.5,axis.args=list(cex.axis=1.5))
    box()
    axis(1, at = seq(0,1650,by=500),labels= c(0,500,1000,1500),cex.axis = 1.5)
    axis(2, at = seq(1650,0,by=-500),labels= c(0,500,1000,1500),cex.axis = 1.5)
    title(paste("Sub-",subs[i]," - high gamma",sep=""),cex.main = 2)
    dev.off()
    
    # collate results
    dBallpowerltg[,,i] <- dBpowerltg
    dBallthetapowerltg[,,i] <- dBthetapowerltg
    dBallalphapowerltg[,,i] <- dBalphapowerltg
    dBallbetapowerltg[,,i] <- dBbetapowerltg
    dBallgammapowerltg[,,i] <- dBgammapowerltg
    dBallhighgammapowerltg[,,i] <- dBhighgammapowerltg
    rm(dBpowerltg,dBthetapowerltg,dBalphapowerltg,dBbetapowerltg,dBgammapowerltg,dBhighgammapowerltg)
    
    # get coefficients
    
    # to find number of electrodes, find the number of coefficients in an
    # arbitrary window, subtract one (for the intercept), and divide by number
    # of frequencies
    tmp <- dBpowermodels[[1]]
    nelecs <- (dim(as.matrix(coef(tmp)))[1]-1)/60
    
    # initialise
    dBpowercoefficients <- array(0,dim=c(nelecs,ntpts,60))
    dBthetapowercoefficients <- array(0,dim=c(nelecs,ntpts,11))
    dBalphapowercoefficients <- array(0,dim=c(nelecs,ntpts,7))
    dBbetapowercoefficients <- array(0,dim=c(nelecs,ntpts,13))
    dBgammapowercoefficients <- array(0,dim=c(nelecs,ntpts,10))
    dBhighgammapowercoefficients<- array(0,dim=c(nelecs,ntpts,19))

    # fill in matrices
    for (j in c(1:ntpts)){
      
      coefs <- as.matrix(coef(dBpowermodels[[j]]))
      # drop the intercept
      coefs <- coefs[-1]
      # reshape to be electrodes x frequencies
      coefs <- t(matrix(coefs,nrow=60,ncol=nelecs))
      dBpowercoefficients[,j,] <- coefs
      
      coefs <- as.matrix(coef(dBthetapowermodels[[j]]))
      # drop the intercept
      coefs <- coefs[-1]
      # reshape to be electrodes x frequencies
      coefs <- t(matrix(coefs,nrow=11,ncol=nelecs))
      dBthetapowercoefficients[,j,] <- coefs
      
      coefs <- as.matrix(coef(dBalphapowermodels[[j]]))
      # drop the intercept
      coefs <- coefs[-1]
      # reshape to be electrodes x frequencies
      coefs <- t(matrix(coefs,nrow=7,ncol=nelecs))
      dBalphapowercoefficients[,j,] <- coefs
      
      coefs <- as.matrix(coef(dBbetapowermodels[[j]]))
      # drop the intercept
      coefs <- coefs[-1]
      # reshape to be electrodes x frequencies
      coefs <- t(matrix(coefs,nrow=13,ncol=nelecs))
      dBbetapowercoefficients[,j,] <- coefs
      
      coefs <- as.matrix(coef(dBgammapowermodels[[j]]))
      # drop the intercept
      coefs <- coefs[-1]
      # reshape to be electrodes x frequencies
      coefs <- t(matrix(coefs,nrow=10,ncol=nelecs))
      dBgammapowercoefficients[,j,] <- coefs
      
      coefs <- as.matrix(coef(dBhighgammapowermodels[[j]]))
      # drop the intercept
      coefs <- coefs[-1]
      # reshape to be electrodes x frequencies
      coefs <- t(matrix(coefs,nrow=19,ncol=nelecs))
      dBhighgammapowercoefficients[,j,] <- coefs
      
    }
    
    # save 
    save(dBpowercoefficients,file=paste(work,"/dB_power_coefficients.Rdata",sep=""))
    save(dBthetapowercoefficients,file=paste(work,"/dB_theta_power_coefficients.Rdata",sep=""))
    save(dBalphapowercoefficients,file=paste(work,"/dB_alpha_power_coefficients.Rdata",sep=""))
    save(dBbetapowercoefficients,file=paste(work,"/dB_beta_power_coefficients.Rdata",sep=""))
    save(dBgammapowercoefficients,file=paste(work,"/dB_gamma_power_coefficients.Rdata",sep=""))
    save(dBhighgammapowercoefficients,file=paste(work,"/dB_high_gamma_power_coefficients.Rdata",sep=""))
    
    # plot selection over time for each electrode - whole-spectrum power only
    
    for (j in c(1:nelecs)){
      if (fileformat=="png"){png(paste(fpath,"/elec_",j,"_dB_power_coefs.png",sep=""))
      }else if (fileformat=="svg") {svg(paste(fpath,"/elec_",j,"_dB_power_coefs.svg",sep=""))}
      plot.coefs(dBpowercoefficients[j,,])
      dev.off()
    }
    
    # plot selection over time for all electrodes in a single plot. Calculate the mean over electrodes (all we really care about is selection, not magnitude or direction)
    
    if (fileformat=="png"){png(paste(fpath,"/elecs_all_dB_power_coefs.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(fpath,"/elecs_all_dB_power_coefs.svg",sep=""))}
    plot.coefs(apply(dBpowercoefficients,c(2:3),mean))
    dev.off()
    
    rm(dBpowercoefficients,dBthetapowercoefficients,dBalphapowercoefficients,dBbetapowercoefficients,dBgammapowercoefficients,dBhighgammapowercoefficients,dBpowermodels,dBthetapowermodels,dBalphapowermodels,dBbetapowermodels,dBgammapowermodels,dBhighgammapowermodels)
    
    # phase 

    # load all non-normalised power data
    
    load(paste(path,"/sub-",subs[i],"/phase_results.Rdata",sep=""))
    phaseacc <- phase[[1]]
    phaseltg <- phase[[2]]
    phasemodels <- phase[[3]]
    load(paste(path,"/sub-",subs[i],"/theta_phase_results.Rdata",sep=""))
    thetaphaseacc <- thetaphase[[1]]
    thetaphaseltg <- thetaphase[[2]]
    thetaphasemodels <- thetaphase[[3]]
    load(paste(path,"/sub-",subs[i],"/alpha_phase_results.Rdata",sep=""))
    alphaphaseacc <- alphaphase[[1]]
    alphaphaseltg <- alphaphase[[2]]
    alphaphasemodels <- alphaphase[[3]]
    load(paste(path,"/sub-",subs[i],"/beta_phase_results.Rdata",sep=""))
    betaphaseacc <- betaphase[[1]]
    betaphaseltg <- betaphase[[2]]
    betaphasemodels <- betaphase[[3]]
    load(paste(path,"/sub-",subs[i],"/gamma_phase_results.Rdata",sep=""))
    gammaphaseacc <- gammaphase[[1]]
    gammaphaseltg <- gammaphase[[2]]
    gammaphasemodels <- gammaphase[[3]]
    load(paste(path,"/sub-",subs[i],"/high_gamma_phase_results.Rdata",sep=""))
    highgammaphaseacc <- highgammaphase[[1]]
    highgammaphaseltg <- highgammaphase[[2]]
    highgammaphasemodels <- highgammaphase[[3]]
    rm(phase,thetaphase,alphaphase,betaphase,gammaphase,highgammaphase)
    
    # plot all frequencies against voltage
    if (fileformat=="png"){png(paste(fpath,"/phase_freq_vs_volt.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(fpath,"/phase_freq_vs_volt.svg",sep=""))}
    plot.cis(voltageacc, colour="#996035")
    plot.cis(phaseacc, colour="#0066FF", newflag = F)
    title(paste("Sub-",subs[i],sep=""), cex.main = 2)
    # add stars to indicate difference from chance
    stars(d1=voltageacc,colour1="#996035",yval=1)
    stars(d1=phaseacc,colour1="#0066FF",yval=0.975)
    # add stars to indicate difference from each other
    stars(d1=voltageacc,d2=phaseacc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
    #save
    dev.off()
    
    # collate results
    groupphaseacc[i,] <- colMeans(phaseacc)
    groupvoltageacc[i,] <- colMeans(voltageacc)
    rm(phaseacc)
    
    # plot all bands
    if (fileformat=="png"){png(paste(fpath,"/phase_bands.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(fpath,"/phase_bands.svg",sep=""))}
    plot.cis(thetaphaseacc, colour="#A50026")
    stars(d1=thetaphaseacc,colour1="#A50026",yval=1)
    plot.cis(alphaphaseacc, colour="#F46D43", newflag = F)
    stars(d1=alphaphaseacc,colour1="#F46D43",yval=0.975)
    plot.cis(betaphaseacc, colour="#FDCC3F", newflag = F)
    stars(d1=betaphaseacc,colour1="#FDCC3F",yval=0.95)
    plot.cis(gammaphaseacc, colour="#66BD63", newflag = F)
    stars(d1=gammaphaseacc,colour1="#66BD63",yval=0.925)
    plot.cis(highgammaphaseacc, colour="#006837", newflag = F)
    stars(d1=highgammaphaseacc,colour1="#006837",yval=0.9)
    title(paste("Sub-",subs[i],sep=""), cex.main = 2)
    dev.off()
    
    # collate results
    groupthetaphaseacc[i,] <- colMeans(thetaphaseacc)
    groupalphaphaseacc[i,] <- colMeans(alphaphaseacc)
    groupbetaphaseacc[i,] <- colMeans(betaphaseacc)
    groupgammaphaseacc[i,] <- colMeans(gammaphaseacc)
    grouphighgammaphaseacc[i,] <- colMeans(highgammaphaseacc)
    rm(thetaphaseacc,alphaphaseacc,betaphaseacc,gammaphaseacc,highgammaphaseacc)
    
    # voltage

    # plot local temporal generalisation matrix. Reversing the index means that the axes are oriented correctly
    if (fileformat=="png"){png(paste(fpath,"/voltage_LTG.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(fpath,"/voltage_LTG",sep=""))}
    image.plot(seq(0,1650,by=10),seq(0,1650,by=10),voltageltg[,ntpts:1],zlim=c(0.4,0.8),xaxt="n",yaxt="n",xlab="Test time (ms) →",ylab="← Train time (ms)",cex.lab = 1.5,axis.args=list(cex.axis=1.5))
    box()
    axis(1, at = seq(0,1650,by=500),labels= c(0,500,1000,1500),cex.axis = 1.5)
    axis(2, at = seq(1650,0,by=-500),labels= c(0,500,1000,1500),cex.axis = 1.5)
    title(paste("Sub-",subs[i]," - voltage",sep=""),cex.main = 2)
    dev.off()

    # collate results
    allvoltageltg[,,i] <- voltageltg
    rm(voltageacc,voltageltg)
    
    # get coefficients
    
    # to find number of electrodes, find the number of coefficients in an
    # arbitrary window, subtract one (for the intercept), and divide by number
    # of timepoints
    tmp <- voltagemodels[[1]]
    nelecs <- (dim(as.matrix(coef(tmp)))[1]-1)/50
    
    # initialise
    voltagecoefficients <- array(0,dim=c(nelecs,ntpts,50))
    
    # fill in matrices
    for (j in c(1:ntpts)){
      
      coefs <- as.matrix(coef(voltagemodels[[j]]))
      # drop the intercept
      coefs <- coefs[-1]
      # reshape to be electrodes x frequencies
      coefs <- t(matrix(coefs,nrow=50,ncol=nelecs))
      voltagecoefficients[,j,] <- coefs
      
    }
    
    # save 
    save(voltagecoefficients,file=paste(work,"/voltage_coefficients.Rdata",sep=""))
    rm(voltagecoefficients,voltagemodels)
  }

  # EOG

  for(i in c(1:length(eogsubs))){

    # decibel-normalised power

    load(paste(path,"/sub-",eogsubs[i],"/dB_power_results.Rdata",sep=""))
    dBpoweracc <- dBpower[[1]]
    load(paste(path,"/sub-",eogsubs[i],"/dB_theta_power_results.Rdata",sep=""))
    dBthetapoweracc <- dBthetapower[[1]]
    load(paste(path,"/sub-",eogsubs[i],"/dB_alpha_power_results.Rdata",sep=""))
    dBalphapoweracc <- dBalphapower[[1]]
    load(paste(path,"/sub-",eogsubs[i],"/dB_beta_power_results.Rdata",sep=""))
    dBbetapoweracc <- dBbetapower[[1]]
    load(paste(path,"/sub-",eogsubs[i],"/dB_gamma_power_results.Rdata",sep=""))
    dBgammapoweracc <- dBgammapower[[1]]
    load(paste(path,"/sub-",eogsubs[i],"/dB_high_gamma_power_results.Rdata",sep=""))
    dBhighgammapoweracc <- dBhighgammapower[[1]]
    rm(dBpower,dBthetapower,dBalphapower,dBbetapower,dBgammapower,dBhighgammapower)
    load(paste(path,"/sub-",eogsubs[i],"/eyes_dB_power_results.Rdata",sep=""))
    load(paste(path,"/sub-",eogsubs[i],"/eyes_dB_theta_power_results.Rdata",sep=""))
    load(paste(path,"/sub-",eogsubs[i],"/eyes_dB_alpha_power_results.Rdata",sep=""))
    load(paste(path,"/sub-",eogsubs[i],"/eyes_dB_beta_power_results.Rdata",sep=""))
    load(paste(path,"/sub-",eogsubs[i],"/eyes_dB_gamma_power_results.Rdata",sep=""))
    load(paste(path,"/sub-",eogsubs[i],"/eyes_dB_high_gamma_power_results.Rdata",sep=""))

    if (fileformat=="png"){png(paste(figurespath,"/sub-",eogsubs[i],"/dB_power_eyes.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(figurespath,"/sub-",eogsubs[i],"/dB_power_eyes.svg",sep=""))}
    plot.cis(dBpowereyes, colour="#9E9E9E")
    plot.cis(dBpoweracc, colour="#0066FF", newflag = F)
    title(paste("Sub-",eogsubs[i],sep=""), cex.main = 2)
    # add stars to indicate difference from chance
    stars(d1=dBpowereyes,colour1="#9E9E9E",yval=1)
    stars(d1=dBpoweracc,colour1="#0066FF",yval=0.975)
    # add stars to indicate difference from each other
    stars(d1=dBpowereyes,d2=dBpoweracc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
    #save
    dev.off()

    if (fileformat=="png"){png(paste(figurespath,"/sub-",eogsubs[i],"/dB_theta_power_eyes.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(figurespath,"/sub-",eogsubs[i],"/dB_theta_power_eyes.svg",sep=""))}
    plot.cis(dBthetapowereyes, colour="#9E9E9E")
    plot.cis(dBthetapoweracc, colour="#A50026", newflag = F)
    title(paste("Sub-",eogsubs[i],sep=""), cex.main = 2)
    # add stars to indicate difference from chance
    stars(d1=dBthetapowereyes,colour1="#9E9E9E",yval=1)
    stars(d1=dBthetapoweracc,colour1="#A50026",yval=0.975)
    # add stars to indicate difference from each other
    stars(d1=dBthetapowereyes,d2=dBthetapoweracc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
    #save
    dev.off()

    if (fileformat=="png"){png(paste(figurespath,"/sub-",eogsubs[i],"/dB_alpha_power_eyes.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(figurespath,"/sub-",eogsubs[i],"/dB_alpha_power_eyes.svg",sep=""))}
    plot.cis(dBalphapowereyes, colour="#9E9E9E")
    plot.cis(dBalphapoweracc, colour="#F46D43", newflag = F)
    title(paste("Sub-",eogsubs[i],sep=""), cex.main = 2)
    # add stars to indicate difference from chance
    stars(d1=dBalphapowereyes,colour1="#9E9E9E",yval=1)
    stars(d1=dBalphapoweracc,colour1="#F46D43",yval=0.975)
    # add stars to indicate difference from each other
    stars(d1=dBalphapowereyes,d2=dBalphapoweracc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
    #save
    dev.off()

    if (fileformat=="png"){png(paste(figurespath,"/sub-",eogsubs[i],"/dB_beta_power_eyes.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(figurespath,"/sub-",eogsubs[i],"/dB_beta_power_eyes.svg",sep=""))}
    plot.cis(dBbetapowereyes, colour="#9E9E9E")
    plot.cis(dBbetapoweracc, colour="#FDCC3F", newflag = F)
    title(paste("Sub-",eogsubs[i],sep=""), cex.main = 2)
    # add stars to indicate difference from chance
    stars(d1=dBbetapowereyes,colour1="#9E9E9E",yval=1)
    stars(d1=dBbetapoweracc,colour1="#FDCC3F",yval=0.975)
    # add stars to indicate difference from each other
    stars(d1=dBbetapowereyes,d2=dBbetapoweracc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
    #save
    dev.off()

    if (fileformat=="png"){png(paste(figurespath,"/sub-",eogsubs[i],"/dB_gamma_power_eyes.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(figurespath,"/sub-",eogsubs[i],"/dB_gamma_power_eyes.svg",sep=""))}
    plot.cis(dBgammapowereyes, colour="#9E9E9E")
    plot.cis(dBgammapoweracc, colour="#66BD63", newflag = F)
    title(paste("Sub-",eogsubs[i],sep=""), cex.main = 2)
    # add stars to indicate difference from chance
    stars(d1=dBgammapowereyes,colour1="#9E9E9E",yval=1)
    stars(d1=dBgammapoweracc,colour1="#66BD63",yval=0.975)
    # add stars to indicate difference from each other
    stars(d1=dBgammapowereyes,d2=dBgammapoweracc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
    #save
    dev.off()

    if (fileformat=="png"){png(paste(figurespath,"/sub-",eogsubs[i],"/dB_high_gamma_power_eyes.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(figurespath,"/sub-",eogsubs[i],"/dB_high_gamma_power_eyes.svg",sep=""))}
    plot.cis(dBhighgammapowereyes, colour="#9E9E9E")
    plot.cis(dBhighgammapoweracc, colour="#006837", newflag = F)
    title(paste("Sub-",eogsubs[i],sep=""), cex.main = 2)
    # add stars to indicate difference from chance
    stars(d1=dBhighgammapowereyes,colour1="#9E9E9E",yval=1)
    stars(d1=dBhighgammapoweracc,colour1="#006837",yval=0.975)
    # add stars to indicate difference from each other
    stars(d1=dBhighgammapowereyes,d2=dBhighgammapoweracc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
    #save
    dev.off()

    # collate results
    dBgrouppowereyes[i,] <- colMeans(dBpowereyes)
    dBgroupthetapowereyes[i,] <- colMeans(dBthetapowereyes)
    dBgroupalphapowereyes[i,] <- colMeans(dBalphapowereyes)
    dBgroupbetapowereyes[i,] <- colMeans(dBbetapowereyes)
    dBgroupgammapowereyes[i,] <- colMeans(dBgammapowereyes)
    dBgrouphighgammapowereyes[i,] <- colMeans(dBhighgammapowereyes)
    rm(dBpowereyes,dBthetapowereyes,dBalphapowereyes,dBbetapowereyes,dBgammapowereyes,dBhighgammapowereyes)

    # 3. PHASE

    load(paste(path,"/sub-",eogsubs[i],"/phase_results.Rdata",sep=""))
    phaseacc <- phase[[1]]
    load(paste(path,"/sub-",eogsubs[i],"/theta_phase_results.Rdata",sep=""))
    thetaphaseacc <- thetaphase[[1]]
    load(paste(path,"/sub-",eogsubs[i],"/alpha_phase_results.Rdata",sep=""))
    alphaphaseacc <- alphaphase[[1]]
    load(paste(path,"/sub-",eogsubs[i],"/beta_phase_results.Rdata",sep=""))
    betaphaseacc <- betaphase[[1]]
    load(paste(path,"/sub-",eogsubs[i],"/gamma_phase_results.Rdata",sep=""))
    gammaphaseacc <- gammaphase[[1]]
    load(paste(path,"/sub-",eogsubs[i],"/high_gamma_phase_results.Rdata",sep=""))
    highgammaphaseacc <- highgammaphase[[1]]
    load(paste(path,"/sub-",eogsubs[i],"/eyes_phase_results.Rdata",sep=""))
    load(paste(path,"/sub-",eogsubs[i],"/eyes_theta_phase_results.Rdata",sep=""))
    load(paste(path,"/sub-",eogsubs[i],"/eyes_alpha_phase_results.Rdata",sep=""))
    load(paste(path,"/sub-",eogsubs[i],"/eyes_beta_phase_results.Rdata",sep=""))
    load(paste(path,"/sub-",eogsubs[i],"/eyes_gamma_phase_results.Rdata",sep=""))
    load(paste(path,"/sub-",eogsubs[i],"/eyes_high_gamma_phase_results.Rdata",sep=""))

    if (fileformat=="png"){png(paste(figurespath,"/sub-",eogsubs[i],"/phase_eyes.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(figurespath,"/sub-",eogsubs[i],"/phase_eyes.svg",sep=""))}
    plot.cis(phaseeyes, colour="#9E9E9E")
    plot.cis(phaseacc, colour="#0066FF", newflag = F)
    title(paste("Sub-",eogsubs[i],sep=""), cex.main = 2)
    # add stars to indicate difference from chance
    stars(d1=phaseeyes,colour1="#9E9E9E",yval=1)
    stars(d1=phaseacc,colour1="#0066FF",yval=0.975)
    # add stars to indicate difference from each other
    stars(d1=phaseeyes,d2=phaseacc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
    #save
    dev.off()
    
    if (fileformat=="png"){png(paste(figurespath,"/sub-",eogsubs[i],"/theta_phase_eyes.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(figurespath,"/sub-",eogsubs[i],"/theta_phase_eyes.svg",sep=""))}
    plot.cis(thetaphaseeyes, colour="#9E9E9E")
    plot.cis(thetaphaseacc, colour="#A50026", newflag = F)
    title(paste("Sub-",eogsubs[i],sep=""), cex.main = 2)
    # add stars to indicate difference from chance
    stars(d1=thetaphaseeyes,colour1="#9E9E9E",yval=1)
    stars(d1=thetaphaseacc,colour1="#A50026",yval=0.975)
    # add stars to indicate difference from each other
    stars(d1=thetaphaseeyes,d2=thetaphaseacc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
    #save
    dev.off()
    
    if (fileformat=="png"){png(paste(figurespath,"/sub-",eogsubs[i],"/alpha_phase_eyes.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(figurespath,"/sub-",eogsubs[i],"/alpha_phase_eyes.svg",sep=""))}
    plot.cis(alphaphaseeyes, colour="#9E9E9E")
    plot.cis(alphaphaseacc, colour="#F46D43", newflag = F)
    title(paste("Sub-",eogsubs[i],sep=""), cex.main = 2)
    # add stars to indicate difference from chance
    stars(d1=alphaphaseeyes,colour1="#9E9E9E",yval=1)
    stars(d1=alphaphaseacc,colour1="#F46D43",yval=0.975)
    # add stars to indicate difference from each other
    stars(d1=alphaphaseeyes,d2=alphaphaseacc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
    #save
    dev.off()
    
    if (fileformat=="png"){png(paste(figurespath,"/sub-",eogsubs[i],"/beta_phase_eyes.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(figurespath,"/sub-",eogsubs[i],"/beta_phase_eyes.svg",sep=""))}
    plot.cis(betaphaseeyes, colour="#9E9E9E")
    plot.cis(betaphaseacc, colour="#FDCC3F", newflag = F)
    title(paste("Sub-",eogsubs[i],sep=""), cex.main = 2)
    # add stars to indicate difference from chance
    stars(d1=betaphaseeyes,colour1="#9E9E9E",yval=1)
    stars(d1=betaphaseacc,colour1="#FDCC3F",yval=0.975)
    # add stars to indicate difference from each other
    stars(d1=betaphaseeyes,d2=betaphaseacc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
    #save
    dev.off()
    
    if (fileformat=="png"){png(paste(figurespath,"/sub-",eogsubs[i],"/gamma_phase_eyes.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(figurespath,"/sub-",eogsubs[i],"/gamma_phase_eyes.svg",sep=""))}
    plot.cis(gammaphaseeyes, colour="#9E9E9E")
    plot.cis(gammaphaseacc, colour="#66BD63", newflag = F)
    title(paste("Sub-",eogsubs[i],sep=""), cex.main = 2)
    # add stars to indicate difference from chance
    stars(d1=gammaphaseeyes,colour1="#9E9E9E",yval=1)
    stars(d1=gammaphaseacc,colour1="#66BD63",yval=0.975)
    # add stars to indicate difference from each other
    stars(d1=gammaphaseeyes,d2=gammaphaseacc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
    #save
    dev.off()
    
    if (fileformat=="png"){png(paste(figurespath,"/sub-",eogsubs[i],"/high_gamma_phase_eyes.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(figurespath,"/sub-",eogsubs[i],"/high_gamma_phase_eyes.svg",sep=""))}
    plot.cis(highgammaphaseeyes, colour="#9E9E9E")
    plot.cis(highgammaphaseacc, colour="#006837", newflag = F)
    title(paste("Sub-",eogsubs[i],sep=""), cex.main = 2)
    # add stars to indicate difference from chance
    stars(d1=highgammaphaseeyes,colour1="#9E9E9E",yval=1)
    stars(d1=highgammaphaseacc,colour1="#006837",yval=0.975)
    # add stars to indicate difference from each other
    stars(d1=highgammaphaseeyes,d2=highgammaphaseacc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
    #save
    dev.off()
    
    # collate results
    groupphaseeyes[i,] <- colMeans(phaseeyes)
    groupthetaphaseeyes[i,] <- colMeans(thetaphaseeyes)
    groupalphaphaseeyes[i,] <- colMeans(alphaphaseeyes)
    groupbetaphaseeyes[i,] <- colMeans(betaphaseeyes)
    groupgammaphaseeyes[i,] <- colMeans(gammaphaseeyes)
    grouphighgammaphaseeyes[i,] <- colMeans(highgammaphaseeyes)
    rm(dBpowereyes,dBthetapowereyes,dBalphapowereyes,dBbetapowereyes,dBgammapowereyes,dBhighgammapowereyes)
    

    # 4. VOLTAGE

    load(paste(path,"/sub-",eogsubs[i],"/voltage_results.Rdata",sep=""))
    voltageacc <- voltage[[1]]
    load(paste(path,"/sub-",eogsubs[i],"/eyes_voltage_results.Rdata",sep=""))

    if (fileformat=="png"){png(paste(figurespath,"/sub-",eogsubs[i],"/voltage_eyes.png",sep=""))
    }else if (fileformat=="svg") {svg(paste(figurespath,"/sub-",eogsubs[i],"/voltage_eyes.svg",sep=""))}
    plot.cis(voltageeyes, colour="#9E9E9E")
    plot.cis(voltageacc, colour="#996035", newflag = F)
    title(paste("Sub-",eogsubs[i],sep=""), cex.main = 2)
    # add stars to indicate difference from chance
    stars(d1=voltageeyes,colour1="#9E9E9E",yval=1)
    stars(d1=voltageacc,colour1="#996035",yval=0.975)
    # add stars to indicate difference from each other
    stars(d1=voltageeyes,d2=voltageacc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
    #save
    dev.off()

    # collate results
    groupvoltageeyes[i,] <- colMeans(voltageeyes)
    rm(voltageeyes)

  }

  # save group-level results. Local temporal generalisation matrices are divided by the number of participants
  dBgrouppowerddp <- list(dBgrouppowerrawacc,dBgrouppowerfilteredacc,dBgrouppowerchannelsrejectedacc,dBgrouppowerrereferencedacc,dBgrouppowertrialsrejectedacc)
  save(dBgrouppowerddp,file=paste(grouppath,"dB_group_power_ddp.Rdata",sep=""))
  groupphaseddp <- list(groupphaserawacc,groupphasefilteredacc,groupphasechannelsrejectedacc,groupphaserereferencedacc,groupphasetrialsrejectedacc)
  save(groupphaseddp,file=paste(grouppath,"group_phase_ddp.Rdata",sep=""))
  groupvoltageddp <- list(groupvoltagerawacc,groupvoltagefilteredacc,groupvoltagechannelsrejectedacc,groupvoltagerereferencedacc,groupvoltagetrialsrejectedacc)
  save(groupvoltageddp,file=paste(grouppath,"group_voltage_ddp.Rdata",sep=""))
  save(dBgrouppoweracc,file=paste(grouppath,"dB_group_power_results.Rdata",sep=""))
  save(dBgroupthetapoweracc,file=paste(grouppath,"dB_group_theta_power_results.Rdata",sep=""))
  save(dBgroupalphapoweracc,file=paste(grouppath,"dB_group_alpha_power_results.Rdata",sep=""))
  save(dBgroupbetapoweracc,file=paste(grouppath,"dB_group_beta_power_results.Rdata",sep=""))
  save(dBgroupgammapoweracc,file=paste(grouppath,"dB_group_gamma_power_results.Rdata",sep=""))
  save(dBgrouphighgammapoweracc,file=paste(grouppath,"dB_group_high_gamma_power_results.Rdata",sep=""))
  save(groupphaseacc,file=paste(grouppath,"group_phase_results.Rdata",sep=""))
  save(groupthetaphaseacc,file=paste(grouppath,"group_theta_phase_results.Rdata",sep=""))
  save(groupalphaphaseacc,file=paste(grouppath,"group_alpha_phase_results.Rdata",sep=""))
  save(groupbetaphaseacc,file=paste(grouppath,"group_beta_phase_results.Rdata",sep=""))
  save(groupgammaphaseacc,file=paste(grouppath,"group_gamma_phase_results.Rdata",sep=""))
  save(grouphighgammaphaseacc,file=paste(grouppath,"group_high_gamma_phase_results.Rdata",sep=""))
  save(groupvoltageacc,file=paste(grouppath,"group_voltage_results.Rdata",sep=""))
  save(dBgrouppowereyes,file=paste(grouppath,"dB_eyes_group_power_results.Rdata",sep=""))
  save(dBgroupthetapowereyes,file=paste(grouppath,"dB_eyes_group_theta_power_results.Rdata",sep=""))
  save(dBgroupalphapowereyes,file=paste(grouppath,"dB_eyes_group_alpha_power_results.Rdata",sep=""))
  save(dBgroupbetapowereyes,file=paste(grouppath,"dB_eyes_group_beta_power_results.Rdata",sep=""))
  save(dBgroupgammapowereyes,file=paste(grouppath,"dB_eyes_group_gamma_power_results.Rdata",sep=""))
  save(dBgrouphighgammapowereyes,file=paste(grouppath,"dB_eyes_group_high_gamma_power_results.Rdata",sep=""))
  dBgrouppowerltg <- apply(dBallpowerltg,c(1:2),mean)
  save(dBgrouppowerltg,file=paste(grouppath,"dB_group_power_LTG.Rdata",sep=""))
  dBgroupthetapowerltg <- apply(dBallthetapowerltg,c(1:2),mean)
  save(dBgroupthetapowerltg,file=paste(grouppath,"dB_group_thetapower_LTG.Rdata",sep=""))
  dBgroupalphapowerltg <- apply(dBallalphapowerltg,c(1:2),mean)
  save(dBgroupalphapowerltg,file=paste(grouppath,"dB_group_alpha_power_LTG.Rdata",sep=""))
  dBgroupbetapowerltg <- apply(dBallbetapowerltg,c(1:2),mean)
  save(dBgroupbetapowerltg,file=paste(grouppath,"dB_group_beta_power_LTG.Rdata",sep=""))
  dBgroupgammapowerltg <- apply(dBallgammapowerltg,c(1:2),mean)
  save(dBgroupgammapowerltg,file=paste(grouppath,"dB_group_gamma_power_LTG.Rdata",sep=""))
  dBgrouphighgammapowerltg <- apply(dBallhighgammapowerltg,c(1:2),mean)
  save(dBgrouphighgammapowerltg,file=paste(grouppath,"dB_group_high_gamma_power_LTG.Rdata",sep=""))
  groupphaseltg <- apply(allphaseltg,c(1:2),mean)
  save(groupphaseltg,file=paste(grouppath,"group_phase_LTG.Rdata",sep=""))
  groupthetaphaseltg <- apply(allthetaphaseltg,c(1:2),mean)
  save(groupthetaphaseltg,file=paste(grouppath,"group_thetaphase_LTG.Rdata",sep=""))
  groupalphaphaseltg <- apply(allalphaphaseltg,c(1:2),mean)
  save(groupalphaphaseltg,file=paste(grouppath,"group_alpha_phase_LTG.Rdata",sep=""))
  groupbetaphaseltg <- apply(allbetaphaseltg,c(1:2),mean)
  save(groupbetaphaseltg,file=paste(grouppath,"group_beta_phase_LTG.Rdata",sep=""))
  groupgammaphaseltg <- apply(allgammaphaseltg,c(1:2),mean)
  save(groupgammaphaseltg,file=paste(grouppath,"group_gamma_phase_LTG.Rdata",sep=""))
  grouphighgammaphaseltg <- apply(allhighgammaphaseltg,c(1:2),mean)
  save(grouphighgammaphaseltg,file=paste(grouppath,"group_high_gamma_phase_LTG.Rdata",sep=""))
  groupvoltageltg <- apply(allvoltageltg,c(1:2),mean)
  save(groupvoltageltg,file=paste(grouppath,"group_voltage_LTG.Rdata",sep=""))
  
  ###### GROUP-LEVEL ANALYSIS #####
  

  # decode during preprocessing
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_power_DDP.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_power_DDP.svg",sep=""))}
  plot.cis(dBgrouppowerrawacc, colour="#583C88")
  stars(d1=dBgrouppowerrawacc,colour1="#583C88",yval=1)
  plot.cis(dBgrouppowerfilteredacc, colour="#9B4FA5", newflag = F)
  stars(d1=dBgrouppowerfilteredacc,colour1="#9B4FA5",yval=0.975)
  plot.cis(dBgrouppowerchannelsrejectedacc, colour="#C4619E", newflag = F)
  stars(d1=dBgrouppowerchannelsrejectedacc,colour1="#C4619E",yval=0.95)
  plot.cis(dBgrouppowerrereferencedacc, colour="#DE89B0", newflag = F)
  stars(d1=dBgrouppowerrereferencedacc,colour1="#DE89B0",yval=0.925)
  plot.cis(dBgrouppowertrialsrejectedacc, colour="#F7BCCB", newflag = F)
  stars(d1=dBgrouppowertrialsrejectedacc,colour1="#F7BCCB",yval=0.9)
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/phase_DDP.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/phase_DDP.svg",sep=""))}
  plot.cis(groupphaserawacc, colour="#583C88")
  stars(d1=groupphaserawacc,colour1="#583C88",yval=1)
  plot.cis(groupphasefilteredacc, colour="#9B4FA5", newflag = F)
  stars(d1=groupphasefilteredacc,colour1="#9B4FA5",yval=0.975)
  plot.cis(groupphasechannelsrejectedacc, colour="#C4619E", newflag = F)
  stars(d1=groupphasechannelsrejectedacc,colour1="#C4619E",yval=0.95)
  plot.cis(groupphaserereferencedacc, colour="#DE89B0", newflag = F)
  stars(d1=groupphaserereferencedacc,colour1="#DE89B0",yval=0.925)
  plot.cis(groupphasetrialsrejectedacc, colour="#F7BCCB", newflag = F)
  stars(d1=groupphasetrialsrejectedacc,colour1="#F7BCCB",yval=0.9)
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/voltage_DDP.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/voltage_DDP.svg",sep=""))}
  plot.cis(groupvoltagerawacc, colour="#583C88")
  stars(d1=groupvoltagerawacc,colour1="#583C88",yval=1)
  plot.cis(groupvoltagefilteredacc, colour="#9B4FA5", newflag = F)
  stars(d1=groupvoltagefilteredacc,colour1="#9B4FA5",yval=0.975)
  plot.cis(groupvoltagechannelsrejectedacc, colour="#C4619E", newflag = F)
  stars(d1=groupvoltagechannelsrejectedacc,colour1="#C4619E",yval=0.95)
  plot.cis(groupvoltagerereferencedacc, colour="#DE89B0", newflag = F)
  stars(d1=groupvoltagerereferencedacc,colour1="#DE89B0",yval=0.925)
  plot.cis(groupvoltagetrialsrejectedacc, colour="#F7BCCB", newflag = F)
  stars(d1=groupvoltagetrialsrejectedacc,colour1="#F7BCCB",yval=0.9)
  dev.off()
  
  # plot ECoG against EOG - (power, phase, voltage) - paired t-tests

  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_power_eyes.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_power_eyes.svg",sep=""))}
  tmppower <- dBgrouppoweracc[c(1,2,3,6,18),]
  plot.cis(dBgrouppowereyes, colour="#9E9E9E")
  plot.cis(tmppower, colour="#0066FF", newflag = F)
  # add stars to indicate difference from chance
  stars(d1=dBgrouppowereyes,colour1="#9E9E9E",yval=1)
  stars(d1=tmppower,colour1="#0066FF",yval=0.975)
  # add stars to indicate difference from each other
  stars(d1=dBgrouppowereyes,d2=tmppower,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_theta_power_eyes.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_theta_power_eyes.svg",sep=""))}
  tmpthetapower <- dBgroupthetapoweracc[c(1,2,3,6,18),]
  plot.cis(dBgroupthetapowereyes, colour="#9E9E9E")
  plot.cis(tmpthetapower, colour="#A50026", newflag = F)
  # add stars to indicate difference from chance
  stars(d1=dBgroupthetapowereyes,colour1="#9E9E9E",yval=1)
  stars(d1=tmpthetapower,colour1="#A50026",yval=0.975)
  # add stars to indicate difference from each other
  stars(d1=dBgroupthetapowereyes,d2=tmpthetapower,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()

  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_alpha_power_eyes.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_alpha_power_eyes.svg",sep=""))}
  tmpalphapower <- dBgroupalphapoweracc[c(1,2,3,6,18),]
  plot.cis(dBgroupalphapowereyes, colour="#9E9E9E")
  plot.cis(tmpalphapower, colour="#F46D43", newflag = F)
  # add stars to indicate difference from chance
  stars(d1=dBgroupalphapowereyes,colour1="#9E9E9E",yval=1)
  stars(d1=tmpalphapower,colour1="#F46D43",yval=0.975)
  # add stars to indicate difference from each other
  stars(d1=dBgroupalphapowereyes,d2=tmpalphapower,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_beta_power_eyes.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_beta_power_eyes.svg",sep=""))}
  tmpbetapower <- dBgroupbetapoweracc[c(1,2,3,6,18),]
  plot.cis(dBgroupbetapowereyes, colour="#9E9E9E")
  plot.cis(tmpbetapower, colour="#FDCC3F", newflag = F)
  # add stars to indicate difference from chance
  stars(d1=dBgroupbetapowereyes,colour1="#9E9E9E",yval=1)
  stars(d1=tmpbetapower,colour1="#FDCC3F",yval=0.975)
  # add stars to indicate difference from each other
  stars(d1=dBgroupbetapowereyes,d2=tmpbetapower,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()

  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_gamma_power_eyes.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_gamma_power_eyes.svg",sep=""))}
  tmpgammapower <- dBgroupgammapoweracc[c(1,2,3,6,18),]
  plot.cis(dBgroupgammapowereyes, colour="#9E9E9E")
  plot.cis(tmpgammapower, colour="#66BD63", newflag = F)
  # add stars to indicate difference from chance
  stars(d1=dBgroupgammapowereyes,colour1="#9E9E9E",yval=1)
  stars(d1=tmpgammapower,colour1="#66BD63",yval=0.975)
  # add stars to indicate difference from each other
  stars(d1=dBgroupgammapowereyes,d2=tmpgammapower,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()

  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_high_gamma_power_eyes.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_high_gamma_power_eyes.svg",sep=""))}
  tmphighgammapower <- dBgrouphighgammapoweracc[c(1,2,3,6,18),]
  plot.cis(dBgrouphighgammapowereyes, colour="#9E9E9E")
  plot.cis(tmphighgammapower, colour="#006837", newflag = F)
  # add stars to indicate difference from chance
  stars(d1=dBgrouphighgammapowereyes,colour1="#9E9E9E",yval=1)
  stars(d1=tmphighgammapower,colour1="#006837",yval=0.975)
  # add stars to indicate difference from each other
  stars(d1=dBgrouphighgammapowereyes,d2=tmphighgammapower,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()

  if (fileformat=="png"){png(paste(groupfigurespath,"/phase_eyes.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/phase_eyes.svg",sep=""))}
  tmpphase <- groupphaseacc[c(1,2,3,6,18),]
  plot.cis(groupphaseeyes, colour="#9E9E9E")
  plot.cis(tmpphase, colour="#0066FF", newflag = F)
  # add stars to indicate difference from chance
  stars(d1=groupphaseeyes,colour1="#9E9E9E",yval=1)
  stars(d1=tmpphase,colour1="#0066FF",yval=0.975)
  # add stars to indicate difference from each other
  stars(d1=groupphaseeyes,d2=tmpphase,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/theta_phase_eyes.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/theta_phase_eyes.svg",sep=""))}
  tmpthetaphase <- groupthetaphaseacc[c(1,2,3,6,18),]
  plot.cis(groupthetaphaseeyes, colour="#9E9E9E")
  plot.cis(tmpthetaphase, colour="#A50026", newflag = F)
  # add stars to indicate difference from chance
  stars(d1=groupthetaphaseeyes,colour1="#9E9E9E",yval=1)
  stars(d1=tmpthetaphase,colour1="#A50026",yval=0.975)
  # add stars to indicate difference from each other
  stars(d1=groupthetaphaseeyes,d2=tmpthetaphase,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/alpha_phase_eyes.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/alph_phase_eyes.svg",sep=""))}
  tmpalphaphase <- groupalphaphaseacc[c(1,2,3,6,18),]
  plot.cis(groupalphaphaseeyes, colour="#9E9E9E")
  plot.cis(tmpalphaphase, colour="#F46D43", newflag = F)
  # add stars to indicate difference from chance
  stars(d1=groupalphaphaseeyes,colour1="#9E9E9E",yval=1)
  stars(d1=tmpalphaphase,colour1="#F46D43",yval=0.975)
  # add stars to indicate difference from each other
  stars(d1=groupalphaphaseeyes,d2=tmpalphaphase,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/beta_phase_eyes.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/beta_phase_eyes.svg",sep=""))}
  tmpbetaphase <- groupbetaphaseacc[c(1,2,3,6,18),]
  plot.cis(groupbetaphaseeyes, colour="#9E9E9E")
  plot.cis(tmpbetaphase, colour="#FDCC3F", newflag = F)
  # add stars to indicate difference from chance
  stars(d1=groupbetaphaseeyes,colour1="#9E9E9E",yval=1)
  stars(d1=tmpbetaphase,colour1="#FDCC3F",yval=0.975)
  # add stars to indicate difference from each other
  stars(d1=groupbetaphaseeyes,d2=tmpbetaphase,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/gamma_phase_eyes.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/gamma_phase_eyes.svg",sep=""))}
  tmpgammaphase <- groupgammaphaseacc[c(1,2,3,6,18),]
  plot.cis(groupgammaphaseeyes, colour="#9E9E9E")
  plot.cis(tmpgammaphase, colour="#66BD63", newflag = F)
  # add stars to indicate difference from chance
  stars(d1=groupgammaphaseeyes,colour1="#9E9E9E",yval=1)
  stars(d1=tmpgammaphase,colour1="#66BD63",yval=0.975)
  # add stars to indicate difference from each other
  stars(d1=groupgammaphaseeyes,d2=tmpgammaphase,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/high_gamma_phase_eyes.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/high_gamma_phase_eyes.svg",sep=""))}
  tmphighgammaphase <- grouphighgammaphaseacc[c(1,2,3,6,18),]
  plot.cis(grouphighgammaphaseeyes, colour="#9E9E9E")
  plot.cis(tmphighgammaphase, colour="#006837", newflag = F)
  # add stars to indicate difference from chance
  stars(d1=grouphighgammaphaseeyes,colour1="#9E9E9E",yval=1)
  stars(d1=tmphighgammaphase,colour1="#006837",yval=0.975)
  # add stars to indicate difference from each other
  stars(d1=grouphighgammaphaseeyes,d2=tmphighgammaphase,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/voltage_eyes.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/voltage_eyes.svg",sep=""))}
  tmpvoltage <- groupvoltageacc[c(1,2,3,6,18),]
  plot.cis(groupvoltageeyes, colour="#9E9E9E")
  plot.cis(tmpvoltage, colour="#996035", newflag = F)
  # add stars to indicate difference from chance
  stars(d1=groupvoltageeyes,colour1="#9E9E9E",yval=1)
  stars(d1=tmpvoltage,colour1="#996035",yval=0.975)
  # add stars to indicate difference from each other
  stars(d1=groupvoltageeyes,d2=tmpvoltage,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()

  
  # Are we justified in combining results from old and new participants and from
  # both hemispheres?
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_power_old_LH_freq_vs_volt.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_old_LH_freq_vs_volt.svg",sep=""))}
  tmpvoltage <- groupvoltageacc[c(1,2,3,4,5,7,9,10),]
  tmpdBpower <- dBgrouppoweracc[c(1,2,3,4,5,7,9,10),]
  plot.cis(tmpvoltage, colour="#996035")
  plot.cis(tmpdBpower, colour="#0066FF", newflag = F)
  # add stars to indicate difference from chance
  stars(d1=tmpvoltage,colour1="#996035",yval=1)
  stars(d1=tmpdBpower,colour1="#0066FF",yval=0.975)
  # add stars to indicate difference from each other
  stars(d1=tmpvoltage,d2=tmpdBpower,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/phase_old_LH_freq_vs_volt.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/phase_old_LH_freq_vs_volt.svg",sep=""))}
  tmpvoltage <- groupvoltageacc[c(1,2,3,4,5,7,9,10),]
  tmpphase <- groupphaseacc[c(1,2,3,4,5,7,9,10),]
  plot.cis(tmpvoltage, colour="#996035")
  plot.cis(tmpphase, colour="#0066FF", newflag = F)
  # add stars to indicate difference from chance
  stars(d1=tmpvoltage,colour1="#996035",yval=1)
  stars(d1=tmpphase,colour1="#0066FF",yval=0.975)
  # add stars to indicate difference from each other
  stars(d1=tmpvoltage,d2=tmpphase,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()
  
  # plot (good) new left hemisphere participants

  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_power_new_LH_freq_vs_volt.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_new_LH_freq_vs_volt.svg",sep=""))}
  tmpvoltage <- groupvoltageacc[c(11,13,14,15,16,17,18),]
  tmpdBpower <- dBgrouppoweracc[c(11,13,14,15,16,17,18),]
  plot.cis(tmpvoltage, colour="#996035")
  plot.cis(tmpdBpower, colour="#0066FF", newflag = F)
  # add stars to indicate difference from chance
  stars(d1=tmpvoltage,colour1="#996035",yval=1)
  stars(d1=tmpdBpower,colour1="#0066FF",yval=0.975)
  # add stars to indicate difference from each other
  stars(d1=tmpvoltage,d2=tmpdBpower,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/phase_new_LH_freq_vs_volt.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/phase_new_LH_freq_vs_volt.svg",sep=""))}
  tmpvoltage <- groupvoltageacc[c(11,13,14,15,16,17,18),]
  tmpphase <- groupphaseacc[c(11,13,14,15,16,17,18),]
  plot.cis(tmpvoltage, colour="#996035")
  plot.cis(tmpphase, colour="#0066FF", newflag = F)
  # add stars to indicate difference from chance
  stars(d1=tmpvoltage,colour1="#996035",yval=1)
  stars(d1=tmpphase,colour1="#0066FF",yval=0.975)
  # add stars to indicate difference from each other
  stars(d1=tmpvoltage,d2=tmpphase,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()
  
  # plot right hemisphere participants
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_power_RH_freq_vs_volt.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_RH_freq_vs_volt.svg",sep=""))}
  tmpvoltage <- groupvoltageacc[c(6,8,12),]
  tmpdBpower <- dBgrouppoweracc[c(6,8,12),]
  plot.cis(tmpvoltage, colour="#996035")
  plot.cis(tmpdBpower, colour="#0066FF", newflag = F)
  # add stars to indicate difference from chance
  stars(d1=tmpvoltage,colour1="#996035",yval=1)
  stars(d1=tmpdBpower,colour1="#0066FF",yval=0.975)
  # add stars to indicate difference from each other
  stars(d1=tmpvoltage,d2=tmpdBpower,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/phase_RH_freq_vs_volt.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/phase_RH_freq_vs_volt.svg",sep=""))}
  tmpvoltage <- groupvoltageacc[c(6,8,12),]
  tmpphase <- groupphaseacc[c(6,8,12),]
  plot.cis(tmpvoltage, colour="#996035")
  plot.cis(tmpphase, colour="#0066FF", newflag = F)
  # add stars to indicate difference from chance
  stars(d1=tmpvoltage,colour1="#996035",yval=1)
  stars(d1=tmpphase,colour1="#0066FF",yval=0.975)
  # add stars to indicate difference from each other
  stars(d1=tmpvoltage,d2=tmpphase,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()
  
  # Because of unbalanced numbers (particularly for the right hemisphere), it
  # is not useful to compare samples with statistics. Nevertheless, we observe
  # the same pattern in all samples - power and voltage perform comparably and 
  # phase performs significantly worse than voltage.
  
  ###################################
  
  # 2. MAIN ANALYSIS 
  
  # plot all frequencies against voltage (dBpower, phase) - paired t-tests and one-sample t-tests against chance
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_power_freq_vs_volt.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_power_freq_vs_volt.svg",sep=""))}
  plot.cis(groupvoltageacc, colour="#996035")
  plot.cis(dBgrouppoweracc, colour="#0066FF", newflag = F)
  # add stars to indicate difference from chance
  stars(d1=groupvoltageacc,colour1="#996035",yval=1)
  stars(d1=dBgrouppoweracc,colour1="#0066FF",yval=0.975)
  # add stars to indicate difference from each other
  stars(d1=groupvoltageacc,d2=dBgrouppoweracc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/phase_freq_vs_volt.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/phase_freq_vs_volt.svg",sep=""))}
  plot.cis(groupvoltageacc, colour="#996035")
  plot.cis(groupphaseacc, colour="#0066FF", newflag = F)
  # add stars to indicate difference from chance
  stars(d1=groupvoltageacc,colour1="#996035",yval=1)
  stars(d1=groupphaseacc,colour1="#0066FF",yval=0.975)
  # add stars to indicate difference from each other
  stars(d1=groupvoltageacc,d2=groupphaseacc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()
  
  # plot all bands (dBpower, phase) - one-sample t-tests against chance (0.5)
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_power_bands.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_power_bands",sep=""))}
  plot.cis(dBgroupthetapoweracc, colour="#A50026")
  stars(d1=dBgroupthetapoweracc,colour1="#A50026",yval=1)
  plot.cis(dBgroupalphapoweracc, colour="#F46D43", newflag = F)
  stars(d1=dBgroupalphapoweracc,colour1="#F46D43",yval=0.975)
  plot.cis(dBgroupbetapoweracc, colour="#FDCC3F", newflag = F)
  stars(d1=dBgroupbetapoweracc,colour1="#FDCC3F",yval=0.95)
  plot.cis(dBgroupgammapoweracc, colour="#66BD63", newflag = F)
  stars(d1=dBgroupgammapoweracc,colour1="#66BD63",yval=0.925)
  plot.cis(dBgrouphighgammapoweracc, colour="#006837", newflag = F)
  stars(d1=dBgrouphighgammapoweracc,colour1="#006837",yval=0.9)
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/phase_bands.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/phase_bands",sep=""))}
  plot.cis(groupthetaphaseacc, colour="#A50026")
  stars(d1=groupthetaphaseacc,colour1="#A50026",yval=1)
  plot.cis(groupalphaphaseacc, colour="#F46D43", newflag = F)
  stars(d1=groupalphaphaseacc,colour1="#F46D43",yval=0.975)
  plot.cis(groupbetaphaseacc, colour="#FDCC3F", newflag = F)
  stars(d1=groupbetaphaseacc,colour1="#FDCC3F",yval=0.95)
  plot.cis(groupgammaphaseacc, colour="#66BD63", newflag = F)
  stars(d1=groupgammaphaseacc,colour1="#66BD63",yval=0.925)
  plot.cis(grouphighgammaphaseacc, colour="#006837", newflag = F)
  stars(d1=grouphighgammaphaseacc,colour1="#006837",yval=0.9)
  dev.off()
  
  # plot power in each band against all frequencies - paired t-tests and one-sample t-tests against chance
  
  # theta
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_theta_power_compare.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_theta_power_compare.svg",sep=""))}
  plot.cis(dBgrouppoweracc, colour="#0066FF")
  plot.cis(dBgroupthetapoweracc, colour="#A50026", newflag = F)
  stars(d1=dBgrouppoweracc,colour1="#0066FF",yval=1)
  stars(d1=dBgroupthetapoweracc,colour1="#A50026",yval=0.975)
  stars(d1=dBgrouppoweracc,d2=dBgroupthetapoweracc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()

  # alpha

  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_alpha_power_compare.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_alpha_power_compare.svg",sep=""))}
  plot.cis(dBgrouppoweracc, colour="#0066FF")
  plot.cis(dBgroupalphapoweracc, colour="#F46D43", newflag = F)
  stars(d1=dBgrouppoweracc,colour1="#0066FF",yval=1)
  stars(d1=dBgroupalphapoweracc,colour1="#F46D43",yval=0.975)
  stars(d1=dBgrouppoweracc,d2=dBgroupalphapoweracc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()
  
  # beta
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_beta_power_compare.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_beta_power_compare.svg",sep=""))}
  plot.cis(dBgrouppoweracc, colour="#0066FF")
  plot.cis(dBgroupbetapoweracc, colour="#FDCC3F", newflag = F)
  stars(d1=dBgrouppoweracc,colour1="#0066FF",yval=1)
  stars(d1=dBgroupbetapoweracc,colour1="#FDCC3F",yval=0.975)
  stars(d1=dBgrouppoweracc,d2=dBgroupbetapoweracc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()
  
  # gamma
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_gamma_power_compare.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_gamma_power_compare.svg",sep=""))}
  plot.cis(dBgrouppoweracc, colour="#0066FF")
  plot.cis(dBgroupgammapoweracc, colour="#66BD63", newflag = F)
  stars(d1=dBgrouppoweracc,colour1="#0066FF",yval=1)
  stars(d1=dBgroupgammapoweracc,colour1="#66BD63",yval=0.975)
  stars(d1=dBgrouppoweracc,d2=dBgroupgammapoweracc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()
  
  # high gamma
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_high_gamma_power_compare.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_high_gamma_power_compare.svg",sep=""))}
  plot.cis(dBgrouppoweracc, colour="#0066FF")
  plot.cis(dBgrouphighgammapoweracc, colour="#006837", newflag = F)
  stars(d1=dBgrouppoweracc,colour1="#0066FF",yval=1)
  stars(d1=dBgrouphighgammapoweracc,colour1="#006837",yval=0.975)
  stars(d1=dBgrouppoweracc,d2=dBgrouphighgammapoweracc,colour1="#000000",colour2="#000000",yval = 0.95,pair=T)
  dev.off()


  # plot local temporal generalisation
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_power_LTG.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_power_LTG.svg",sep=""))}
  image.plot(seq(0,1650,by=10),seq(0,1650,by=10),dBgrouppowerltg[,ntpts:1],zlim=c(0.4,0.8),xaxt="n",yaxt="n",xlab="Test time (ms) →",ylab="← Train time (ms)",cex.lab = 1.5,axis.args=list(cex.axis=1.5))
  box()
  axis(1, at = seq(0,1650,by=500),labels= c(0,500,1000,1500),cex.axis = 1.5)
  axis(2, at = seq(1650,0,by=-500),labels= c(0,500,1000,1500),cex.axis = 1.5)
  dev.off()

  
  if (fileformat=="png"){png(paste(groupfigurespath,"/theta_dBpower_LTG.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/theta_dBpower_LTG.svg",sep=""))}
  image.plot(seq(0,1650,by=10),seq(0,1650,by=10),dBgroupthetapowerltg[,ntpts:1],zlim=c(0.4,0.8),xaxt="n",yaxt="n",xlab="Test time (ms) →",ylab="← Train time (ms)",cex.lab = 1.5,axis.args=list(cex.axis=1.5))
  box()
  axis(1, at = seq(0,1650,by=500),labels= c(0,500,1000,1500),cex.axis = 1.5)
  axis(2, at = seq(1650,0,by=-500),labels= c(0,500,1000,1500),cex.axis = 1.5)
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/alpha_dBpower_LTG.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/alpha_dBpower_LTG.svg",sep=""))}
  image.plot(seq(0,1650,by=10),seq(0,1650,by=10),dBgroupalphapowerltg[,ntpts:1],zlim=c(0.4,0.8),xaxt="n",yaxt="n",xlab="Test time (ms) →",ylab="← Train time (ms)",cex.lab = 1.5,axis.args=list(cex.axis=1.5))
  box()
  axis(1, at = seq(0,1650,by=500),labels= c(0,500,1000,1500),cex.axis = 1.5)
  axis(2, at = seq(1650,0,by=-500),labels= c(0,500,1000,1500),cex.axis = 1.5)
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/beta_dBpower_LTG.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/beta_dBpower_LTG.svg",sep=""))}
  image.plot(seq(0,1650,by=10),seq(0,1650,by=10),dBgroupbetapowerltg[,ntpts:1],zlim=c(0.4,0.8),xaxt="n",yaxt="n",xlab="Test time (ms) →",ylab="← Train time (ms)",cex.lab = 1.5,axis.args=list(cex.axis=1.5))
  box()
  axis(1, at = seq(0,1650,by=500),labels= c(0,500,1000,1500),cex.axis = 1.5)
  axis(2, at = seq(1650,0,by=-500),labels= c(0,500,1000,1500),cex.axis = 1.5)
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/gamma_dBpower_LTG.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/gamma_dBpower_LTG.svg",sep=""))}
  image.plot(seq(0,1650,by=10),seq(0,1650,by=10),dBgroupgammapowerltg[,ntpts:1],zlim=c(0.4,0.8),xaxt="n",yaxt="n",xlab="Test time (ms) →",ylab="← Train time (ms)",cex.lab = 1.5,axis.args=list(cex.axis=1.5))
  box()
  axis(1, at = seq(0,1650,by=500),labels= c(0,500,1000,1500),cex.axis = 1.5)
  axis(2, at = seq(1650,0,by=-500),labels= c(0,500,1000,1500),cex.axis = 1.5)
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/high_gamma_dBpower_LTG.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/high_gamma_dBpower_LTG.svg",sep=""))}
  image.plot(seq(0,1650,by=10),seq(0,1650,by=10),dBgrouphighgammapowerltg[,ntpts:1],zlim=c(0.4,0.8),xaxt="n",yaxt="n",xlab="Test time (ms) →",ylab="← Train time (ms)",cex.lab = 1.5,axis.args=list(cex.axis=1.5))
  box()
  axis(1, at = seq(0,1650,by=500),labels= c(0,500,1000,1500),cex.axis = 1.5)
  axis(2, at = seq(1650,0,by=-500),labels= c(0,500,1000,1500),cex.axis = 1.5)
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/voltage_LTG.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/voltage_LTG.svg",sep=""))}
  image.plot(seq(0,1650,by=10),seq(0,1650,by=10),groupvoltageltg[,ntpts:1],zlim=c(0.4,0.8),xaxt="n",yaxt="n",xlab="Test time (ms) →",ylab="← Train time (ms)",cex.lab = 1.5,axis.args=list(cex.axis=1.5))
  box()
  axis(1, at = seq(0,1650,by=500),labels= c(0,500,1000,1500),cex.axis = 1.5)
  axis(2, at = seq(1650,0,by=-500),labels= c(0,500,1000,1500),cex.axis = 1.5)
  dev.off()
  
  # plot local temporal generalisation matrix for all bands, averaged together
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_power_average_bands__LTG.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_power_average_bands_LTG.svg",sep=""))}
  tmp <- (dBgroupthetapowerltg + dBgroupalphapowerltg + dBgroupbetapowerltg + dBgroupgammapowerltg + dBgroupgammapowerltg)/5
  image.plot(seq(0,1650,by=10),seq(0,1650,by=10),tmp[,ntpts:1],zlim=c(0.4,0.8),xaxt="n",yaxt="n",xlab="Test time (ms) →",ylab="← Train time (ms)",cex.lab = 1.5,axis.args=list(cex.axis=1.5))
  box()
  axis(1, at = seq(0,1650,by=500),labels= c(0,500,1000,1500),cex.axis = 1.5)
  axis(2, at = seq(1650,0,by=-500),labels= c(0,500,1000,1500),cex.axis = 1.5)
  dev.off()
  
  # statistical threshold for generalisation plots
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_power_overlap_wav.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_power_overlap_wav.svg",sep=""))}
  plot.cluster.compare(dBallpowerltg,dBallpowerltg,colour1="#D8E7FF",colour2="#0066FF")
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_theta_power_overlap_wav.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_theta_power_overlap_wav.svg",sep=""))}
  plot.cluster.compare(dBallthetapowerltg,dBallthetapowerltg,colour1="#FF7E9C",colour2="#A50026")
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_alpha_power_overlap_wav.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_alpha_power_overlap_wav.svg",sep=""))}
  plot.cluster.compare(dBallalphapowerltg,dBallalphapowerltg,colour1="#FCD2C5",colour2="#F46D43")
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_beta_power_overlap_wav.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_beta_power_overlap_wav.svg",sep=""))}
  plot.cluster.compare(dBallbetapowerltg,dBallbetapowerltg,colour1="#FEF0C7",colour2="#EBaF02")
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_gamma_power_overlap_wav.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_gamma_power_overlap_wav.svg",sep=""))}
  plot.cluster.compare(dBallgammapowerltg,dBallgammapowerltg,colour1="#E0F1DF",colour2="#66BD63")
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_high_gamma_power_overlap_wav.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_high_gamma_power_overlap_wav.svg",sep=""))}
  plot.cluster.compare(dBallhighgammapowerltg,dBallhighgammapowerltg,colour1="#81C97F",colour2="#006837")
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/voltage_overlap_wav.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/voltage_overlap_wav.svg",sep=""))}
  plot.cluster.compare(allvoltageltg,allvoltageltg,colour1="#F2E5DB",colour2="#996035")
  dev.off()
  
  
  # widening generalisation window using 10 clusters
  
  # decibel-normalised power
  
  #get cluster indices from the group-level clusters
  clusters <- get.clust(dBgrouppowerltg,k=10,dmeth="cos")
  
  # initialise array of clustered timecourses
  dBgrouppowerclusterbysub <- array(0,dim=c(max(clusters),ntpts,length(subs)))
  
  # apply clustering indices from group-level clustering to individual participants
  for (i in c(1:max(clusters))){
    # if there is only one classifier in a cluster, don't average
    if (length(which(clusters==i))==1){
      dBgrouppowerclusterbysub[i,,] <- dBallpowerltg[clusters==i,,]
    } else {dBgrouppowerclusterbysub[i,,] <- apply(dBallpowerltg[clusters==i,,], c(2,3), mean)
    }
  }
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_power_WGW_10clust.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_power_WGW_10clust.svg",sep=""))}
  plot.clacc(dBgrouppowerltg,clusters,dBgrouppowerclusterbysub,colour1="#0066FF",colour2="#D8E7FF")
  dev.off()
  
  # theta
  
  #get cluster indices from the group-level clusters
  clusters <- get.clust(dBgroupthetapowerltg,k=10,dmeth="cos")
  
  # initialise array of clustered timecourses
  dBgroupthetapowerclusterbysub <- array(0,dim=c(max(clusters),ntpts,length(subs)))
  
  # apply clustering indices from group-level clustering to individual participants
  for (i in c(1:max(clusters))){
    # if there is only one classifier in a cluster, don't average
    if (length(which(clusters==i))==1){
      dBgroupthetapowerclusterbysub[i,,] <- dBallthetapowerltg[clusters==i,,]
    } else {dBgroupthetapowerclusterbysub[i,,] <- apply(dBallthetapowerltg[clusters==i,,], c(2,3), mean)
    }
  }
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_theta_power_WGW_10clust.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_theta_power_WGW_10clust.svg",sep=""))}
  plot.clacc(dBgroupthetapowerltg,clusters,dBgroupthetapowerclusterbysub,colour1="#A50026",colour2="#FF7E9C")
  dev.off()
  
  # alpha
  
  #get cluster indices from the group-level clusters
  clusters <- get.clust(dBgroupalphapowerltg,k=10,dmeth="cos")
  
  # initialise array of clustered timecourses
  dBgroupalphapowerclusterbysub <- array(0,dim=c(max(clusters),ntpts,length(subs)))
  
  # apply clustering indices from group-level clustering to individual participants
  for (i in c(1:max(clusters))){
    # if there is only one classifier in a cluster, don't average
    if (length(which(clusters==i))==1){
      dBgroupalphapowerclusterbysub[i,,] <- dBallalphapowerltg[clusters==i,,]
    } else {dBgroupalphapowerclusterbysub[i,,] <- apply(dBallalphapowerltg[clusters==i,,], c(2,3), mean)
    }
  }
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_alpha_power_WGW_10clust.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_alpha_power_WGW_10clust.svg",sep=""))}
  plot.clacc(dBgroupalphapowerltg,clusters,dBgroupalphapowerclusterbysub,colour1="#F46D43",colour2="#FCD2C5")
  dev.off()
  
  # beta
  
  #get cluster indices from the group-level clusters
  clusters <- get.clust(dBgroupbetapowerltg,k=10,dmeth="cos")
  
  # initialise array of clustered timecourses
  dBgroupbetapowerclusterbysub <- array(0,dim=c(max(clusters),ntpts,length(subs)))
  
  # apply clustering indices from group-level clustering to individual participants
  for (i in c(1:max(clusters))){
    # if there is only one classifier in a cluster, don't average
    if (length(which(clusters==i))==1){
      dBgroupbetapowerclusterbysub[i,,] <- dBallbetapowerltg[clusters==i,,]
    } else {dBgroupbetapowerclusterbysub[i,,] <- apply(dBallbetapowerltg[clusters==i,,], c(2,3), mean)
    }
  }
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_beta_power_WGW_10clust.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_beta_power_WGW_10clust.svg",sep=""))}
  plot.clacc(dBgroupbetapowerltg,clusters,dBgroupbetapowerclusterbysub,colour1="#EBaF02",colour2="#FEF0C7")
  dev.off()
  
  # gamma
  
  #get cluster indices from the group-level clusters
  clusters <- get.clust(dBgroupgammapowerltg,k=10,dmeth="cos")
  
  # initialise array of clustered timecourses
  dBgroupgammapowerclusterbysub <- array(0,dim=c(max(clusters),ntpts,length(subs)))
  
  # apply clustering indices from group-level clustering to individual participants
  for (i in c(1:max(clusters))){
    # if there is only one classifier in a cluster, don't average
    if (length(which(clusters==i))==1){
      dBgroupgammapowerclusterbysub[i,,] <- dBallgammapowerltg[clusters==i,,]
    } else {dBgroupgammapowerclusterbysub[i,,] <- apply(dBallgammapowerltg[clusters==i,,], c(2,3), mean)
    }
  }
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_gamma_power_WGW_10clust.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_gamma_power_WGW_10clust.svg",sep=""))}
  plot.clacc(dBgroupgammapowerltg,clusters,dBgroupgammapowerclusterbysub,colour1="#66BD63",colour2="#E0F1DF")
  dev.off()
  
  # high gamma
  
  #get cluster indices from the group-level clusters
  clusters <- get.clust(dBgrouphighgammapowerltg,k=10,dmeth="cos")
  
  # initialise array of clustered timecourses
  dBgrouphighgammapowerclusterbysub <- array(0,dim=c(max(clusters),ntpts,length(subs)))
  
  # apply clustering indices from group-level clustering to individual participants
  for (i in c(1:max(clusters))){
    # if there is only one classifier in a cluster, don't average
    if (length(which(clusters==i))==1){
      dBgrouphighgammapowerclusterbysub[i,,] <- dBallhighgammapowerltg[clusters==i,,]
    } else {dBgrouphighgammapowerclusterbysub[i,,] <- apply(dBallhighgammapowerltg[clusters==i,,], c(2,3), mean)
    }
  }
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_high_gamma_power_WGW_10clust.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_high_gamma_power_WGW_10clust.svg",sep=""))}
  plot.clacc(dBgrouphighgammapowerltg,clusters,dBgrouphighgammapowerclusterbysub,colour1="#006837",colour2="#81C97F")
  dev.off()

  # voltage
  
  #get cluster indices from the group-level clusters
  clusters <- get.clust(groupvoltageltg,k=10,dmeth="cos")
  
  # initialise array of clustered timecourses
  voltageclusterbysub <- array(0,dim=c(max(clusters),ntpts,length(subs)))
  
  # apply clustering indices from group-level clustering to individual participants
  for (i in c(1:max(clusters))){
    # if there is only one classifier in a cluster, don't average
    if (length(which(clusters==i))==1){
      voltageclusterbysub[i,,] <- allvoltageltg[clusters==i,,]
    } else {voltageclusterbysub[i,,] <- apply(allvoltageltg[clusters==i,,], c(2,3), mean)
    }
  }
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/voltage_WGW_10clust.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/voltage_WGW_10clust.svg",sep=""))}
  plot.clacc(groupvoltageltg,clusters,voltageclusterbysub,colour1="#996035",colour2="#F2E5DB")
  dev.off()
  
  # plot area under the cluster timecourse
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_power_AUC_widen.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_power_AUC_widen.svg",sep=""))}
  plot.auc.widen(x=seq(0,1650,by=10),ymat=dBgrouppowerltg,colour="#0066FF")
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_theta_power_AUC_widen.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_theta_power_AUC_widen.svg",sep=""))}
  plot.auc.widen(x=seq(0,1650,by=10),ymat=dBgroupthetapowerltg,colour="#A50026")
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_alpha_power_AUC_widen.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_alpha_power_AUC_widen.svg",sep=""))}
  plot.auc.widen(x=seq(0,1650,by=10),ymat=dBgroupalphapowerltg,colour="#F46D43")
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_beta_power_AUC_widen.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_beta_power_AUC_widen.svg",sep=""))}
  plot.auc.widen(x=seq(0,1650,by=10),ymat=dBgroupbetapowerltg,colour="#EBaF02")
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_gamma_power_AUC_widen.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_gamma_power_AUC_widen.svg",sep=""))}
  plot.auc.widen(x=seq(0,1650,by=10),ymat=dBgroupgammapowerltg,colour="#66BD63")
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/dB_high_gamma_power_AUC_widen.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/dB_high_gamma_power_AUC_widen.svg",sep=""))}
  plot.auc.widen(x=seq(0,1650,by=10),ymat=dBgrouphighgammapowerltg,colour="#006837")
  dev.off()
  
  if (fileformat=="png"){png(paste(groupfigurespath,"/voltage_AUC_widen.png",sep=""))
  }else if (fileformat=="svg") {svg(paste(groupfigurespath,"/voltage_AUC_widen.svg",sep=""))}
  plot.auc.widen(x=seq(0,1650,by=10),ymat=groupvoltageltg,colour="#996035")
  dev.off()
  

  # collate coefficients

  work <- gsub("derivatives/modelfit","work",path)
  # load mni coordinates
  setwd("Z:/Saskia/ECoG_LASSO/scripts")
  mni <- read.csv("mni_coordinates.csv",header=T)
  # initialise. As in the mni coordinates file, we stack the electrodes from
  # each participant
  dBgrouppowercoefficients <- array(0,dim=c(dim(mni)[1],ntpts,60))
  dBgroupthetapowercoefficients <- array(0,dim=c(dim(mni)[1],ntpts,11))
  dBgroupalphapowercoefficients <- array(0,dim=c(dim(mni)[1],ntpts,7))
  dBgroupbetapowercoefficients <- array(0,dim=c(dim(mni)[1],ntpts,13))
  dBgroupgammapowercoefficients <- array(0,dim=c(dim(mni)[1],ntpts,10))
  dBgrouphighgammapowercoefficients <- array(0,dim=c(dim(mni)[1],ntpts,19))
  groupvoltagecoefficients <- array(0,dim=c(dim(mni)[1],ntpts,50))

  counter <- 0

  for (i in c(1:length(subs))){
    # load coefficients
    load(paste(work,"/sub-",subs[i],"/dB_power_coefficients.Rdata",sep=""))
    load(paste(work,"/sub-",subs[i],"/dB_theta_power_coefficients.Rdata",sep=""))
    load(paste(work,"/sub-",subs[i],"/dB_alpha_power_coefficients.Rdata",sep=""))
    load(paste(work,"/sub-",subs[i],"/dB_beta_power_coefficients.Rdata",sep=""))
    load(paste(work,"/sub-",subs[i],"/dB_gamma_power_coefficients.Rdata",sep=""))
    load(paste(work,"/sub-",subs[i],"/dB_high_gamma_power_coefficients.Rdata",sep=""))
    load(paste(work,"/sub-",subs[i],"/voltage_coefficients.Rdata",sep=""))
    # find number of electrodes
    nelecs <- dim(voltagecoefficients)[1]
    # add each participants' data to the group
    dBgrouppowercoefficients[c((counter+1):(counter+nelecs)),,] <- dBpowercoefficients
    dBgroupthetapowercoefficients[c((counter+1):(counter+nelecs)),,] <- dBthetapowercoefficients
    dBgroupalphapowercoefficients[c((counter+1):(counter+nelecs)),,] <- dBalphapowercoefficients
    dBgroupbetapowercoefficients[c((counter+1):(counter+nelecs)),,] <- dBbetapowercoefficients
    dBgroupgammapowercoefficients[c((counter+1):(counter+nelecs)),,] <- dBgammapowercoefficients
    dBgrouphighgammapowercoefficients[c((counter+1):(counter+nelecs)),,] <- dBhighgammapowercoefficients
    groupvoltagecoefficients[c((counter+1):(counter+nelecs)),,] <- voltagecoefficients
    # update row counter
    counter = counter + nelecs
  }

  # save as Rdata
  save(dBgrouppowercoefficients,file=paste(coefwork,"dB_group_power_coefficients.Rdata",sep=""))
  save(dBgroupthetapowercoefficients,file=paste(coefwork,"dB_group_power_coefficients.Rdata",sep=""))
  save(dBgroupalphapowercoefficients,file=paste(coefwork,"dB_group_alpha_power_coefficients.Rdata",sep=""))
  save(dBgroupbetapowercoefficients,file=paste(coefwork,"dB_group_beta_power_coefficients.Rdata",sep=""))
  save(dBgroupgammapowercoefficients,file=paste(coefwork,"dB_group_gamma_power_coefficients.Rdata",sep=""))
  save(dBgrouphighgammapowercoefficients,file=paste(coefwork,"dB_group_high_gamma_power_coefficients.Rdata",sep=""))
  save(groupvoltagecoefficients,file=paste(coefwork,"group_voltage_coefficients.Rdata",sep=""))
  
  # save as .csv
  # calculate mean over frequencies (or, in the case of voltage, timepoints within window) before saving
  
  tmp <- dBgrouppowercoefficients
  tmp[tmp==0] <- NA
  coefs <- apply(tmp,c(1:2),mean,na.rm=TRUE)
  coefs[(is.na(coefs))] <- 0
  write.table(coefs,file=paste(coefwork,"/dB_power_coefficients.csv",sep=""),row.names=F,col.names=F)
  
  tmp <- dBgroupthetapowercoefficients
  tmp[tmp==0] <- NA
  coefs <- apply(tmp,c(1:2),mean,na.rm=TRUE)
  coefs[(is.na(coefs))] <- 0
  write.table(coefs,file=paste(coefwork,"/dB_theta_power_coefficients.csv",sep=""),row.names=F,col.names=F)
  
  tmp <- dBgroupalphapowercoefficients
  tmp[tmp==0] <- NA
  coefs <- apply(tmp,c(1:2),mean,na.rm=TRUE)
  coefs[(is.na(coefs))] <- 0
  write.table(coefs,file=paste(coefwork,"/dB_alpha_power_coefficients.csv",sep=""),row.names=F,col.names=F)
  
  tmp <- dBgroupbetapowercoefficients
  tmp[tmp==0] <- NA
  coefs <- apply(tmp,c(1:2),mean,na.rm=TRUE)
  coefs[(is.na(coefs))] <- 0
  write.table(coefs,file=paste(coefwork,"/dB_beta_power_coefficients.csv",sep=""),row.names=F,col.names=F)
  
  tmp <- dBgroupgammapowercoefficients
  tmp[tmp==0] <- NA
  coefs <- apply(tmp,c(1:2),mean,na.rm=TRUE)
  coefs[(is.na(coefs))] <- 0
  write.table(coefs,file=paste(coefwork,"/dB_gamma_power_coefficients.csv",sep=""),row.names=F,col.names=F)
  
  tmp <- dBgrouphighgammapowercoefficients
  tmp[tmp==0] <- NA
  coefs <- apply(tmp,c(1:2),mean,na.rm=TRUE)
  coefs[(is.na(coefs))] <- 0
  write.table(coefs,file=paste(coefwork,"/dB_high_gamma_power_coefficients.csv",sep=""),row.names=F,col.names=F)
  
  tmp <- groupvoltagecoefficients
  tmp[tmp==0] <- NA
  coefs <- apply(tmp,c(1:2),mean,na.rm=TRUE)
  coefs[(is.na(coefs))] <- 0
  write.table(coefs,file=paste(coefwork,"/voltage_coefficients.csv",sep=""),row.names=F,col.names=F)

}
