#!/usr/bin/env Rscript 
####################
## Collect arguments
args <- commandArgs(TRUE)
## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL
rm(argsL)
## Display help when no all arguments passed or help needed
if("--help" %in% args | is.null(args$arg1)) {
  cat("
      The R Script fit_models.r
      Mandatory arguments:
      --arg1=string         - path to folder containing data .csvs at multiple preprocessing stages
      --help                - print this text
      Example:
      ./sub_decode_during_preprocessing.r --arg1=$dirp/derivatives/decode_during_preprocessing/sub-$q/ \n\n")
  q(save="no")
}

######################################

.libPaths("/group/mlr-lab/Saskia/toolboxes/rlibs")

library(glmnet)
library(doMC)

cpus_to_use <- parallel::detectCores()
cpus_to_use
registerDoMC(cores=cpus_to_use)

get.win <- dget("/group/mlr-lab/Saskia/ECoG_LASSO/scripts/get_win.r")
get.freq.win <- dget("/group/mlr-lab/Saskia/ECoG_LASSO/scripts/get_freq_win.r")
fit.models.during.preprocessing <- dget("/group/mlr-lab/Saskia/ECoG_LASSO/scripts/fit_models_during_preprocessing.r")
decode.during.preprocessing <- dget("/group/mlr-lab/Saskia/ECoG_LASSO/scripts/decode_during_preprocessing.r")

path = paste(args$arg1, sep = "")
decode.during.preprocessing(path)


